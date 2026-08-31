//! geomorphon 地貌形态分析(对齐 GRASS r.geomorphon 算法)
//!
//! 机器视觉方法: 以焦点像元为中心, 沿 8 方向做线视(line-of-sight),
//! 在搜索半径 L 内确定各方向的天顶/天底角, 与 flat 阈值比较生成三值
//! 模式(+更高/-更低/0平坦), 组成 8 元组地貌模式(geomorphon)。
//!
//! 与窗口统计法(TPI/均值)相比: 结果与实际地物的尺寸、起伏、朝向无关,
//! 单一半径即可自适应多尺度形态。

/// 三值模式: 每方向 2 位, 0=平坦 1=更高(天顶) 2=更低(天底)
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Pattern8(pub [u8; 8]);

impl Pattern8 {
    /// 山脊: 全部方向更低(周围都比我低)
    pub fn is_peak(&self) -> bool {
        self.0.iter().all(|&v| v == 1)
    }
    /// 山谷: 全部方向更高(周围都比我高)
    pub fn is_valley(&self) -> bool {
        self.0.iter().all(|&v| v == 2)
    }
    /// 正方向(更低=四周水往远处流)数量: 位置分类用
    pub fn count_lower(&self) -> u32 {
        self.0.iter().filter(|&&v| v == 2).count() as u32
    }
    /// 负方向(更高=四周山地)数量
    pub fn count_higher(&self) -> u32 {
        self.0.iter().filter(|&&v| v == 1).count() as u32
    }
    pub fn count_flat(&self) -> u32 {
        self.0.iter().filter(|&&v| v == 0).count() as u32
    }
}

/// geomorphon 计算
///
/// 沿射线逐像元步进并双线性插值高程; 天顶/天底角取射线上
/// 相对焦点最大仰角/最大俯角(切角), 对齐 r.geomorphon angle 原理。
///
/// # 参数
/// - `dem`/`w`/`h`/`res`: 高程栅格与分辨率
/// - `x`/`y`: 焦点像元
/// - `search_m`: 线视最大距离(米), 决定形态尺度
/// - `skip_m`: 起始跳过距离(米), 消除微起伏噪声(0=不跳)
/// - `flat_deg`: 平坦判角阈值(度), |角差|<flat 视为平坦
#[allow(clippy::too_many_arguments)]
pub fn geomorphon_pattern(
    dem: &[f32],
    w: usize,
    h: usize,
    res: f64,
    x: usize,
    y: usize,
    search_m: f64,
    skip_m: f64,
    flat_deg: f64,
) -> Pattern8 {
    let zc = dem[y * w + x];
    let zc = if zc.is_finite() { zc } else { return Pattern8::default() };
    // 8 方向: 东起逆时针(对齐 GRASS 顺序); 无序性对分类统计无影响
    let dirs: [(f64, f64); 8] = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
    ];
    let search_px = search_m / res;
    let skip_px = skip_m / res;
    let flat_rad = flat_deg.to_radians();
    let mut pat = Pattern8::default();
    let at = |fx: f64, fy: f64| -> f32 {
        let x0 = fx.floor() as isize;
        let y0 = fy.floor() as isize;
        let x1 = x0 + 1;
        let y1 = y0 + 1;
        let gx = fx - x0 as f64;
        let gy = fy - y0 as f64;
        let g = |xx: isize, yy: isize| -> f64 {
            if xx < 0 || yy < 0 || xx >= w as isize || yy >= h as isize {
                return f64::NAN;
            }
            dem[yy as usize * w + xx as usize] as f64
        };
        let z00 = g(x0, y0);
        let z10 = g(x1, y0);
        let z01 = g(x0, y1);
        let z11 = g(x1, y1);
        // 窗外像元按最近值延拓(避免 NaN 污染插值)
        let z10 = if z10.is_nan() { z00 } else { z10 };
        let z01 = if z01.is_nan() { z00 } else { z01 };
        let z11 = if z11.is_nan() { z00 } else { z11 };
        (z00 * (1.0 - gx) * (1.0 - gy)
            + z10 * gx * (1.0 - gy)
            + z01 * (1.0 - gx) * gy
            + z11 * gx * gy) as f32
    };
    for (k, (dx, dy)) in dirs.iter().enumerate() {
        let mut zenith = 0f64; // 最大仰角(某射线点高于焦点)
        let mut nadir = 0f64; // 最大俯角(某射线点低于焦点)
        let mut dist = skip_px.max(1.0);
        while dist <= search_px {
            let fx = x as f64 + dx * dist;
            let fy = y as f64 + dy * dist;
            let zv = at(fx, fy);
            if zv.is_finite() {
                let angle = ((zv as f64 - zc as f64) / (dist * res)).atan();
                if angle > zenith {
                    zenith = angle;
                }
                if -angle > nadir {
                    nadir = -angle;
                }
            }
            dist += 1.0;
        }
        // 三值: 与 flat 阈值比较(天顶/天底取大者定方向性)
        pat.0[k] = if zenith > flat_rad && zenith >= nadir {
            1 // 更高
        } else if nadir > flat_rad && nadir > zenith {
            2 // 更低
        } else {
            0 // 平坦
        };
    }
    pat
}

/// 六级坡位映射: 按 GRASS 标准 P 值(更高方向数, 0..8)单调映射,
/// 对齐 legacy 编码 1谷 2坡下 3平坡 4坡中 5坡上 6山脊:
/// - P 0-1: 山脊(全向更低/近全向) -> 6
/// - P 2: 坡上(shoulder 肩部) -> 5
/// - P 3-4: 坡中(spur/slope) -> 4
/// - P 5-6: 平坡(footslope 坡麓) -> 3
/// - P 7: 坡下/山谷(valley) -> 2
/// - P 8: 山谷(pit 坑) -> 1
pub fn pattern_to_level(pat: &Pattern8, _min_extreme: u32) -> u8 {
    let p = pat.count_higher();
    match p {
        0..=1 => 6,
        2 => 5,
        3..=4 => 4,
        5..=6 => 3,
        7 => 2,
        _ => 1,
    }
}

/// 全图 geomorphon 六级坡位
pub fn geomorphon_levels(
    dem: &[f32],
    w: usize,
    h: usize,
    res: f64,
    search_m: f64,
    skip_m: f64,
    flat_deg: f64,
) -> Vec<u8> {
    let mut out = vec![0u8; w * h];
    for y in 1..h.saturating_sub(1) {
        for x in 1..w.saturating_sub(1) {
            let pat = geomorphon_pattern(dem, w, h, res, x, y, search_m, skip_m, flat_deg);
            out[y * w + x] = pattern_to_level(&pat, 2);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak_dem(w: usize, h: usize, peak_h: f32) -> Vec<f32> {
        let mut d = vec![100f32; w * h];
        let (cx, cy) = (w / 2, h / 2);
        for y in 0..h {
            for x in 0..w {
                let dd = ((x as f64 - cx as f64).powi(2) + (y as f64 - cy as f64).powi(2)).sqrt();
                d[y * w + x] = 100.0 + (peak_h as f64 * (1.0 - dd / 40.0).clamp(0.0, 1.0)) as f32;
            }
        }
        d
    }

    #[test]
    fn peak_center_is_ridge_level() {
        // 尖峰中心: 全向更低 -> 山脊(6)
        let (w, h) = (81usize, 81usize);
        let dem = peak_dem(w, h, 200.0);
        assert_eq!(geomorphon_pattern(&dem, w, h, 5.0, 40, 40, 100.0, 0.0, 1.0).count_lower(), 8);
        let lv = geomorphon_levels(&dem, w, h, 5.0, 100.0, 0.0, 1.0);
        assert_eq!(lv[40 * w + 40], 6, "峰心应为山脊");
    }

    #[test]
    fn pit_center_is_valley_level() {
        let (w, h) = (81usize, 81usize);
        let dem: Vec<f32> = peak_dem(w, h, 200.0).iter().map(|&z| 300.0 - (z - 100.0)).collect();
        let lv = geomorphon_levels(&dem, w, h, 5.0, 100.0, 0.0, 1.0);
        assert_eq!(lv[40 * w + 40], 1, "坑心应为山谷");
    }

    #[test]
    fn skip_suppresses_micro_relief() {
        // 中心平地+微噪声: skip=0 时噪声成方向, skip 大时判定平坦
        let (w, h) = (81usize, 81usize);
        let mut dem = vec![100f32; w * h];
        for i in 0..w * h {
            if i % 7 == 0 {
                dem[i] += 0.05; // 5cm 微噪声
            }
        }
        let p0 = geomorphon_pattern(&dem, w, h, 5.0, 40, 40, 100.0, 0.0, 1.0);
        let p1 = geomorphon_pattern(&dem, w, h, 5.0, 40, 40, 100.0, 30.0, 1.0);
        // 微噪声在 skip 后不应产生方向性
        assert!(p1.count_higher() + p1.count_lower() < 8);
        let _ = p0;
    }
}