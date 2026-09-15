//! 确定性合成 DEM 夹具与断言辅助(无随机成分)。
//!
//! 所有高程函数一律使用像元中心米坐标: `x=(col+0.5)*res`, `y=(row+0.5)*res`,
//! 保证同一场景在不同分辨率下解析一致, 供尺度稳定性测试复用。

use topo_core::input::RasterShape;

/// 平面斜坡: 自南向北线性上升(坡降 0.08), 南缘谷底、北缘山脊。
pub fn planar_slope(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let slope = 0.08;
    let mut z = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = (row as f64 + 0.5) * res;
        for _col in 0..width {
            z.push((900.0 + slope * y) as f32);
        }
    }
    (z, RasterShape { width, height, resolution_m: res })
}

/// V 形谷: 中央谷底、两侧线性上升至边缘山脊; 谷底向北微倾保证排水。
pub fn v_valley(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let mut z = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = (row as f64 + 0.5) * res;
        for col in 0..width {
            let x = (col as f64 + 0.5) * res - cx;
            z.push((800.0 + 0.10 * x.abs() + 0.004 * y) as f32);
        }
    }
    (z, RasterShape { width, height, resolution_m: res })
}

/// 宽坦盆地: 圆形平底加环状山脊围限, 盆底沿 x 微倾提供排水出口。
pub fn broad_basin(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z: Vec<f32> = (0..height)
        .flat_map(|row| {
            (0..width).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r = x.hypot(y);
                let floor = 800.0 + 0.002 * x;
                let rim = ((r - 900.0) / 250.0).clamp(0.0, 1.0);
                (floor + 120.0 * rim * rim) as f32
            })
        })
        .collect();
    (z, RasterShape { width, height, resolution_m: res })
}

/// 闭合洼地: 外环高地包围中央低地且无地表出口(喀斯特/山间闭合盆地原型)。
pub fn closed_pit(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z: Vec<f32> = (0..height)
        .flat_map(|row| {
            (0..width).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r = x.hypot(y);
                let v = if r < 400.0 {
                    810.0 + 2.0 * (r / 400.0) * (r / 400.0)
                } else if r < 800.0 {
                    let t = (r - 400.0) / 400.0;
                    812.0 + 60.0 * t * t
                } else {
                    872.0 + 0.08 * (r - 800.0)
                };
                v as f32
            })
        })
        .collect();
    (z, RasterShape { width, height, resolution_m: res })
}

/// 锥形山丘: 中央峰顶(850 m)加周围平原(700 m)。
pub fn conical_hill(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z: Vec<f32> = (0..height)
        .flat_map(|row| {
            (0..width).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r = x.hypot(y);
                let v = 700.0 + 150.0 * (1.0 - r / 600.0).max(0.0);
                v as f32
            })
        })
        .collect();
    (z, RasterShape { width, height, resolution_m: res })
}

/// 嵌套山脊: 沿 x 方向周期性脊—谷交替(谷底 850 m、脊顶 930 m, 波长 640 m),
/// 整体向北缓倾 0.8% 提供汇流方向。
pub fn nested_ridges(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    const WAVELENGTH_M: f64 = 640.0;
    let mut z = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = (row as f64 + 0.5) * res;
        for col in 0..width {
            let x = (col as f64 + 0.5) * res;
            let ridge = 40.0 * (1.0 - (std::f64::consts::TAU * x / WAVELENGTH_M).cos());
            z.push((850.0 + ridge + 0.008 * y) as f32);
        }
    }
    (z, RasterShape { width, height, resolution_m: res })
}

/// 边缘 NoData: 在锥形山丘场景外圈固定 8 像元范围写入 NaN,
/// 模拟与图幅边界连通的无效区(必须保持无效, 不得内插)。
pub fn with_border_nodata(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let (mut z, shape) = conical_hill(width, height, res);
    let border = 8usize;
    for row in 0..height {
        for col in 0..width {
            if row < border || col < border || row + border >= height || col + border >= width {
                z[row * width + col] = f32::NAN;
            }
        }
    }
    (z, shape)
}
