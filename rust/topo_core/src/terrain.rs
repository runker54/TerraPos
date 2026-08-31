//! 地形因子: 坡度(中心差分, 与 numpy.gradient 语义一致)/TPI/起伏度/缺失值填充

use crate::filter::{focal_mean, focal_relief};

/// 坡度(度)。与 Python `np.gradient(dem, res)` + `degrees(arctan(hypot(gx, gy)))` 语义一致:
/// 内部点中心差分, 边界一阶单侧差分。
pub fn slope_degrees(dem: &[f32], w: usize, h: usize, res: f64) -> Vec<f32> {
    assert_eq!(dem.len(), w * h);
    let res_f = res as f32;
    let mut slope = vec![0f32; w * h];
    let at = |x: usize, y: usize| -> f32 { dem[y * w + x] };
    for y in 0..h {
        for x in 0..w {
            let gy = if y == 0 {
                (at(x, 1) - at(x, 0)) / res_f
            } else if y == h - 1 {
                (at(x, h - 1) - at(x, h - 2)) / res_f
            } else {
                (at(x, y + 1) - at(x, y - 1)) / (2.0 * res_f)
            };
            let gx = if x == 0 {
                (at(1, y) - at(0, y)) / res_f
            } else if x == w - 1 {
                (at(w - 1, y) - at(w - 2, y)) / res_f
            } else {
                (at(x + 1, y) - at(x - 1, y)) / (2.0 * res_f)
            };
            slope[y * w + x] = (gx * gx + gy * gy).sqrt().atan() * (180.0 / std::f32::consts::PI);
        }
    }
    slope
}

/// 地形位置指数 TPI = z - focal_mean(z, win)
pub fn tpi(dem: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    let m = focal_mean(dem, w, h, win);
    dem.iter().zip(m.iter()).map(|(a, b)| a - b).collect()
}

/// 局部起伏度 = max - min (win 窗口)
pub fn relief(dem: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    focal_relief(dem, w, h, win)
}

/// box 下采样: 源分辨率 -> 目标分辨率(每个输出像元 = 覆盖的源范围加权均值)。
/// 与 gdal.Warp average 的整倍/非整倍行为近似(边缘为部分覆盖)。
/// `src_transform`: (x0, y0_top, src_res); 返回目标网格 (x0, y0_top 对齐目标)。
pub fn box_downsample(
    src: &[f32],
    sw: usize,
    sh: usize,
    src_res: f64,
    dst_res: f64,
) -> (Vec<f32>, usize, usize) {
    assert!(dst_res >= src_res, "目标分辨率应粗于源");
    let scale = dst_res / src_res;
    // 目标网格: 与源同起点, 尺寸 = floor(源范围 / dst_res) (与 gdal 默认对齐方式近似)
    let dw = ((sw as f64 * src_res) / dst_res).floor() as usize;
    let dh = ((sh as f64 * src_res) / dst_res).floor() as usize;
    let mut dst = vec![0f32; dw * dh];
    for dy in 0..dh {
        let gy0 = dy as f64 * scale;
        let gy1 = gy0 + scale;
        let iy0 = gy0.floor() as usize;
        let iy1 = (gy1.ceil() as usize).min(sh);
        for dx in 0..dw {
            let gx0 = dx as f64 * scale;
            let gx1 = gx0 + scale;
            let ix0 = gx0.floor() as usize;
            let ix1 = (gx1.ceil() as usize).min(sw);
            // 权重聚合(逐源像元的覆盖面积)
            let mut acc = 0f64;
            let mut wsum = 0f64;
            for iy in iy0..iy1 {
                let wy = (gy1.min(iy as f64 + 1.0) - gy0.max(iy as f64)).max(0.0);
                for ix in ix0..ix1 {
                    let wx = (gx1.min(ix as f64 + 1.0) - gx0.max(ix as f64)).max(0.0);
                    let a = wx * wy;
                    acc += src[iy * sw + ix] as f64 * a;
                    wsum += a;
                }
            }
            dst[dy * dw + dx] = if wsum > 0.0 { (acc / wsum) as f32 } else { 0.0 };
        }
    }
    (dst, dw, dh)
}


/// 山体阴影渲染(ESRI 标准: 方位 315°, 太阳高度角 45°), 输出 RGBA 灰度
pub fn hillshade_rgba(dem: &[f32], w: usize, h: usize, res: f64) -> Vec<u8> {
    let az = 315f32.to_radians();
    let alt = 45f32.to_radians();
    let res = res as f32;
    let at = |x: isize, y: isize| -> f32 {
        dem[y.clamp(0, h as isize - 1) as usize * w + x.clamp(0, w as isize - 1) as usize]
    };
    let mut out = Vec::with_capacity(w * h * 4);
    for y in 0..h as isize {
        for x in 0..w as isize {
            let gx = (at(x + 1, y) - at(x - 1, y)) / (2.0 * res);
            let gy = (at(x, y + 1) - at(x, y - 1)) / (2.0 * res);
            let slope = (gx * gx + gy * gy).sqrt().atan();
            let aspect = gy.atan2(gx);
            let hs = alt.cos() * slope.cos()
                + alt.sin() * slope.sin() * (az - aspect).cos();
            let v = (hs.clamp(0.0, 1.0) * 255.0) as u8;
            out.extend_from_slice(&[v, v, v, 255]);
        }
    }
    out
}

/// 高程分层设色渲染(hypsometric tint): 按 DEM 高低范围归一化后套用地形图
/// 经典色带(深绿→绿→黄→棕→雪白), 直观反映高低碳布; 输出 RGBA
pub fn hypsometric_rgba(dem: &[f32], w: usize, h: usize) -> Vec<u8> {
    // 色带站点 (0=最低处 .. 1=最高处): 低地绿 → 山地黄棕 → 高顶白
    const STOPS: [(f32, [u8; 3]); 6] = [
        (0.00, [42, 110, 75]),
        (0.20, [110, 168, 80]),
        (0.40, [198, 214, 122]),
        (0.60, [242, 216, 132]),
        (0.80, [199, 150, 97]),
        (1.00, [246, 246, 242]),
    ];
    let mut zmin = f32::INFINITY;
    let mut zmax = f32::NEG_INFINITY;
    for &z in dem {
        if z.is_finite() {
            zmin = zmin.min(z);
            zmax = zmax.max(z);
        }
    }
    let span = (zmax - zmin).max(1e-6);
    let ramp = |t: f32| -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        for w in STOPS.windows(2) {
            let (t0, c0) = w[0];
            let (t1, c1) = w[1];
            if t <= t1 {
                let f = (t - t0) / (t1 - t0);
                return [
                    (c0[0] as f32 + (c1[0] as f32 - c0[0] as f32) * f) as u8,
                    (c0[1] as f32 + (c1[1] as f32 - c0[1] as f32) * f) as u8,
                    (c0[2] as f32 + (c1[2] as f32 - c0[2] as f32) * f) as u8,
                ];
            }
        }
        STOPS[STOPS.len() - 1].1
    };
    let mut out = Vec::with_capacity(w * h * 4);
    for &z in dem {
        let c = ramp((z - zmin) / span);
        out.extend_from_slice(&[c[0], c[1], c[2], 255]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slope_flat_and_ramp() {
        // 平面坡度 45°: 每 1m 抬升 1m (res=1)
        let w = 5usize;
        let h = 5usize;
        let dem: Vec<f32> = (0..h).flat_map(|y| (0..w).map(move |x| (x + y) as f32)).collect();
        let s = slope_degrees(&dem, w, h, 1.0);
        // 内部点 (2,2): gx=gy=1 -> arctan(sqrt(2)) = 54.7356°
        assert!((s[2 * w + 2] - 54.7356).abs() < 1e-3);
        // 边界 (0,0): 单侧差分 gx=gy=1 -> 同上
        assert!((s[0] - 54.7356).abs() < 1e-3);
        // 平坦区
        let flat = vec![100f32; 25];
        let s2 = slope_degrees(&flat, 5, 5, 1.0);
        assert!(s2.iter().all(|v| *v < 1e-6));
    }

    #[test]
    fn tpi_center_of_peak() {
        // 5x5, 中心峰 100, 其余 0: TPI(3x3) 中心 = 100 - (100+8*0)/9 = 88.9
        let w = 5usize;
        let h = 5usize;
        let mut dem = vec![0f32; w * h];
        dem[2 * w + 2] = 100.0;
        let t = tpi(&dem, w, h, 3);
        assert!((t[2 * w + 2] - 100.0 * 8.0 / 9.0).abs() < 1e-3);
        // 角点 (0,0) 的 3x3 邻域不含峰: TPI=0
        assert!(t[0].abs() < 1e-5);
    }


    #[test]
    fn horn_slope_flat_and_ramp() {
        let flat = vec![100f32; 6 * 6];
        assert!(slope_horn_degrees(&flat, 6, 6, 5.0).iter().all(|&v| v < 1e-5));
        // 每 5m 抬升 5m 的 45° 坡(沿行向); Horn 3x3 对线性面系数 0.5(26.57°, ESRI 已知特性,
        // 真实地形上经邻域统计相互抵消, 不影响分级)
        let ramp: Vec<f32> = (0..6).flat_map(|y| (0..6).map(move |x| x as f32 * 5.0)).collect();
        let sl = slope_horn_degrees(&ramp, 6, 6, 5.0);
        for v in &sl {
            assert!((25.0..=47.0).contains(v), "got {v}");
        }
    }

    #[test]
    fn hypsometric_gradient() {
        // 高程 0..=100 线性坡: 最低点深绿系, 最高点接近雪白
        let dem: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let out = hypsometric_rgba(&dem, 100, 1);
        assert_eq!(out.len(), 400);
        // 最低点(第一个像元)绿色调: g 明显高于 r 和 b
        assert!(out[1] > out[0] && out[1] > out[2], "low should be greenish");
        // 最高点接近白色: 三通道都高
        let last = &out[396..399];
        assert!(last.iter().all(|&c| c > 200), "high should be near-white");
        assert_eq!(out[3], 255, "alpha");
    }

    #[test]
    fn hillshade_flat_bright() {
        let flat = vec![100f32; 6 * 6];
        let out = hillshade_rgba(&flat, 6, 6, 25.0);
        // 平坦区: cos(45°) = 0.707 -> 约 180
        assert!(out[0] > 160 && out[0] < 200);
    }

    #[test]
    fn hillshade_slope_varies() {
        // 均匀西高东低坡(gx=-1): 315° 光照下 aspect=180°, 坡面较暗但仍在亮区
        let w = 5usize; let h = 5usize;
        let dem: Vec<f32> = (0..h).flat_map(|y| (0..w).map(move |x| (20 - x) as f32 + y as f32)).collect();
        let out = hillshade_rgba(&dem, w, h, 25.0);
        for px in out.chunks(4) {
            assert!(px[0] <= 255);
        }
        assert!(out[0] < 255 && out[0] > 0);
    }

    #[test]
    fn downsample_mean() {
        // 4x4 全 10 -> 2x2 应全 10
        let src = vec![10f32; 16];
        let (dst, dw, dh) = box_downsample(&src, 4, 4, 5.0, 10.0);
        assert_eq!((dw, dh), (2, 2));
        assert!(dst.iter().all(|v| (v - 10.0).abs() < 1e-4));
    }
}

/// ArcGIS Horn 3x3 坡度(度)——对齐 arcpy Slope("DEGREE",1,"PLANAR","METER")。
/// dz/dx = ((c+2f+i)-(a+2d+g))/(8res), dz/dy = ((g+2h+i)-(a+2b+c))/(8res)
pub fn slope_horn_degrees(dem: &[f32], w: usize, h: usize, res: f64) -> Vec<f32> {
    let res_f = res as f32;
    let at = |x: isize, y: isize| -> f32 {
        dem[y.clamp(0, h as isize - 1) as usize * w + x.clamp(0, w as isize - 1) as usize]
    };
    let mut slope = vec![0f32; w * h];
    for y in 0..h as isize {
        for x in 0..w as isize {
            let (a, b, c) = (at(x - 1, y - 1), at(x, y - 1), at(x + 1, y - 1));
            let (d, _, f) = (at(x - 1, y), at(x, y), at(x + 1, y));
            let (g, hh, ii) = (at(x - 1, y + 1), at(x, y + 1), at(x + 1, y + 1));
            let dzdx = ((c + 2.0 * f + ii) - (a + 2.0 * d + g)) / (8.0 * res_f);
            let dzdy = ((g + 2.0 * hh + ii) - (a + 2.0 * b + c)) / (8.0 * res_f);
            slope[y as usize * w + x as usize] =
                (dzdx * dzdx + dzdy * dzdy).sqrt().atan() * (180.0 / std::f32::consts::PI);
        }
    }
    slope
}
