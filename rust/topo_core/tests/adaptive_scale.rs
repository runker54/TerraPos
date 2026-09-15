//! 自适应尺度契约测试: 特征尺度分辨率稳定性、地貌宽度响应与平移不变性。

mod common;

use topo_core::input::RasterShape;
use topo_core::scale::{build_scale_pyramid, scale_bounds_for, usable_scales, BASE_SCALES_M};

const GROWTH: f32 = 0.15;

/// 高斯丘场景(解析确定): z = base + amp * exp(-r^2 / (2 sigma^2))
fn gaussian_hill(size: usize, res: f64, sigma_m: f64, amp: f64, base: f64) -> (Vec<f32>, RasterShape) {
    let cx = size as f64 * res * 0.5;
    let cy = size as f64 * res * 0.5;
    let z: Vec<f32> = (0..size)
        .flat_map(|row| {
            (0..size).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r2 = x * x + y * y;
                (base + amp * (-r2 / (2.0 * sigma_m * sigma_m)).exp()) as f32
            })
        })
        .collect();
    (z, RasterShape { width: size, height: size, resolution_m: res })
}

/// 中心地貌(r < 2sigma)上特征尺度的中位数
fn median_scale_central(z: &[f32], shape: RasterShape, scale_m: &[f32], sigma_m: f64) -> f64 {
    let cx = shape.width as f64 * shape.resolution_m * 0.5;
    let cy = shape.height as f64 * shape.resolution_m * 0.5;
    let mut vals: Vec<f32> = (0..shape.height)
        .flat_map(|row| (0..shape.width).map(move |col| (row, col)))
        .filter(|&(row, col)| {
            let x = (col as f64 + 0.5) * shape.resolution_m - cx;
            let y = (row as f64 + 0.5) * shape.resolution_m - cy;
            x.hypot(y) < 2.0 * sigma_m
        })
        .map(|(row, col)| scale_m[row * shape.width + col])
        .collect();
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vals[vals.len() / 2] as f64
}

/// 几何相同的高斯丘在 5/10/25 m 采样下特征尺度相差不超过一个相邻尺度步
#[test]
fn characteristic_scale_stable_across_resolutions() {
    let mut medians = Vec::new();
    for (size, res) in [(801usize, 5.0f64), (401, 10.0), (161, 25.0)] {
        let (z, shape) = gaussian_hill(size, res, 400.0, 120.0, 800.0);
        let valid = vec![true; z.len()];
        let pyramid = build_scale_pyramid(&z, &valid, shape, GROWTH).unwrap();
        medians.push(median_scale_central(&z, shape, &pyramid.characteristic_scale_m, 400.0));
    }
    let scale_index = |m: f64| -> i64 {
        BASE_SCALES_M
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - m).abs().partial_cmp(&(b.1 - m).abs()).unwrap()
            })
            .map(|(i, _)| i as i64)
            .unwrap()
    };
    let idx: Vec<i64> = medians.iter().map(|&m| scale_index(m)).collect();
    assert!(
        idx.iter().all(|&i| (i - idx[0]).abs() <= 1),
        "特征尺度跨分辨率漂移超过一级: medians={medians:?} idx={idx:?}"
    );
}

/// 宽阔地貌选择比窄小地貌更大的分析尺度
#[test]
fn broad_landform_selects_larger_scale_than_narrow() {
    let (zn, sn) = gaussian_hill(201, 25.0, 150.0, 120.0, 800.0);
    let (zb, sb) = gaussian_hill(201, 25.0, 800.0, 120.0, 800.0);
    let valid = vec![true; zn.len()];
    let pn = build_scale_pyramid(&zn, &valid, sn, GROWTH).unwrap();
    let pb = build_scale_pyramid(&zb, &valid, sb, GROWTH).unwrap();
    let mn = median_scale_central(&zn, sn, &pn.characteristic_scale_m, 150.0);
    let mb = median_scale_central(&zb, sb, &pb.characteristic_scale_m, 800.0);
    assert!(
        mb > mn,
        "宽地貌尺度应大于窄地貌: broad={mb} narrow={mn}"
    );
}

/// 整体高程平移不改变特征尺度(逐像元一致)
#[test]
fn vertical_offset_preserves_selected_scales() {
    let (z, shape) = gaussian_hill(201, 25.0, 400.0, 120.0, 800.0);
    let shifted: Vec<f32> = z.iter().map(|&v| v + 50.0).collect();
    let valid = vec![true; z.len()];
    let p1 = build_scale_pyramid(&z, &valid, shape, GROWTH).unwrap();
    let p2 = build_scale_pyramid(&shifted, &valid, shape, GROWTH).unwrap();
    for (a, b) in p1
        .characteristic_scale_m
        .iter()
        .zip(p2.characteristic_scale_m.iter())
    {
        assert_eq!(a, b, "平移后特征尺度必须逐像元一致");
    }
}

/// 可用尺度裁剪: 小于 5 粗像元或大于短边四分之一的尺度被删除
#[test]
fn usable_scales_respect_extent() {
    let shape = RasterShape { width: 400, height: 300, resolution_m: 25.0 };
    let s = usable_scales(shape);
    // res 25m: 最小尺度 >= 125m 全保留; 短边 300*25=7500m, 上限 1875m -> 4000 删
    assert!(s.contains(&125.0) && s.contains(&1000.0));
    assert!(!s.contains(&4000.0), "4000m 超过短边四分之一应删除");
    // 分辨率 40m: 125m < 5*40=200 -> 125 删
    let s2 = usable_scales(RasterShape { width: 400, height: 300, resolution_m: 40.0 });
    assert!(!s2.contains(&125.0), "小于 5 像元的尺度应删除");
}

/// 有效范围不足以支撑至少三个尺度时明确报错
#[test]
fn rejects_extent_too_small_for_multiscale() {
    let (z, shape) = gaussian_hill(60, 25.0, 200.0, 100.0, 800.0);
    let valid = vec![true; z.len()];
    let err = build_scale_pyramid(&z, &valid, shape, GROWTH).unwrap_err().to_string();
    assert!(err.contains("multiscale"), "{err}");
}

/// 地貌亚类尺度上下限(规格 9.4 精确边界)
#[test]
fn scale_bounds_follow_subclass_ranges() {
    // 低丘(<500m, relief<200)
    assert_eq!(scale_bounds_for(400.0f32, 150.0, 3000.0), (125.0f32, 750.0));
    // 高丘
    assert_eq!(scale_bounds_for(400.0f32, 300.0, 3000.0), (250.0f32, 1500.0));
    // 低山 500-1000m
    assert_eq!(scale_bounds_for(800.0f32, 300.0, 3000.0), (250.0f32, 2000.0));
    // 中山 1000-3500m
    assert_eq!(scale_bounds_for(1500.0f32, 300.0, 3000.0), (500.0f32, 4000.0));
    // 高山 3500-5000m
    assert_eq!(scale_bounds_for(4000.0f32, 300.0, 9000.0), (750.0f32, 6000.0));
    // 极高山 >=5000m, 上限受数据约束
    assert_eq!(scale_bounds_for(5200.0f32, 300.0, 7000.0), (1000.0f32, 7000.0));
}
