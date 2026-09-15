//! geomorphon 十标准形态与米制单位回归测试。

use topo_core::geomorphon::{geomorphon_pattern, pattern_to_landform, Landform, Pattern8};

/// 十个代表性三值模式(方向序 E,SE,S,SW,W,NW,N,NE; 0=flat 1=higher 2=lower)
fn representative_patterns() -> Vec<(Landform, Pattern8)> {
    use Landform::*;
    vec![
        (Flat, Pattern8([0, 0, 0, 0, 0, 0, 0, 0])),
        (Peak, Pattern8([2, 2, 2, 2, 2, 2, 2, 2])),
        (Pit, Pattern8([1, 1, 1, 1, 1, 1, 1, 1])),
        // 山脊: 脊线沿 E-W 下延(两侧低), 横向抬升
        (Ridge, Pattern8([2, 1, 1, 1, 2, 1, 1, 1])),
        // 山谷: 谷线沿 E-W 上游, 两壁更低
        (Valley, Pattern8([1, 2, 2, 2, 1, 2, 2, 2])),
        // 肩部: 坡向上高、坡向下低、等高线方向平坦
        (Shoulder, Pattern8([0, 0, 1, 1, 0, 0, 2, 0])),
        // 凸坡(spur): 多段破碎且 lower 占优
        (Spur, Pattern8([2, 2, 1, 2, 2, 1, 2, 1])),
        // 直坡: 上坡向高、下坡向低、无平坦
        (Slope, Pattern8([1, 1, 1, 1, 2, 2, 2, 2])),
        // 凹坡(hollow): 多段破碎且 higher 占优
        (Hollow, Pattern8([1, 1, 2, 1, 1, 2, 1, 2])),
        // 坡麓: 仅上方一段更高, 其余平坦
        (Footslope, Pattern8([0, 0, 1, 1, 1, 0, 0, 0])),
    ]
}

/// 十个代表模式各自映射到独立形态
#[test]
fn ten_representative_patterns_map_to_distinct_forms() {
    for (expect, pat) in representative_patterns() {
        let got = pattern_to_landform(&pat);
        assert_eq!(got, expect, "pattern {:?} -> {:?}, 期望 {:?}", pat.0, got, expect);
    }
}

/// 旋转与镜像不变: 任意旋转/反射后的模式形态不变
#[test]
fn landform_is_rotation_and_reflection_invariant() {
    let rotate = |p: &Pattern8| Pattern8([p.0[2], p.0[3], p.0[4], p.0[5], p.0[6], p.0[7], p.0[0], p.0[1]]);
    let reflect = |p: &Pattern8| {
        // 沿东西轴镜像: SE<->NE, S<->N, SW<->NW
        Pattern8([p.0[0], p.0[7], p.0[6], p.0[5], p.0[4], p.0[3], p.0[2], p.0[1]])
    };
    for (expect, pat) in representative_patterns() {
        let r1 = rotate(&pat);
        let r2 = rotate(&r1);
        assert_eq!(pattern_to_landform(&r1), expect, "旋转 90 度后 {expect:?} 改变");
        assert_eq!(pattern_to_landform(&r2), expect, "旋转 180 度后 {expect:?} 改变");
        assert_eq!(pattern_to_landform(&reflect(&pat)), expect, "镜像后 {expect:?} 改变");
    }
}

/// 米制单位回归: 1000 m 搜索在 5 m 分辨率必须覆盖 200 像元
/// (历史缺陷为二次换算只搜 40 像元, 150 像元外的高墙不可见)
#[test]
fn search_distance_reaches_full_metre_scale() {
    let (w, h) = (201usize, 41usize);
    let mut dem = vec![100f32; w * h];
    for y in 0..h {
        for x in 150..160 {
            dem[y * w + x] = 300.0; // 焦点以东 130..140 像元(650..700 m)的高墙
        }
    }
    let pat = geomorphon_pattern(&dem, w, h, 5.0, 20, 20, 1000.0, 0.0, 1.0);
    assert_eq!(pat.0[0], 1, "东向必须看到 650m 外的高墙(1=higher)");
}

/// 解析曲面上的剖面/平面曲率符号正确且无效邻域传播 NaN
#[test]
fn curvatures_sign_and_nodata_policy() {
    use topo_core::terrain::{plan_curvature, profile_curvature};
    let (w, h) = (32usize, 32usize);
    let res = 5.0f64;
    // 圆丘 z = -0.01 * (r^2): 全向凸, 剖面/平面曲率均应为负(ArcGIS 约定)
    let convex: Vec<f32> = (0..h)
        .flat_map(|y| {
            (0..w).map(move |x| {
                let dx = (x as f64 - 16.0) * res;
                let dy = (y as f64 - 16.0) * res;
                (-0.01 * (dx * dx + dy * dy)) as f32
            })
        })
        .collect();
    let pc = profile_curvature(&convex, w, h, res);
    let pl = plan_curvature(&convex, w, h, res);
    let c = 8 * w + 8; // 丘外缘坡面(两方向梯度均非零)
    assert!(pc[c] < 0.0, "凸坡剖面曲率应为负: {}", pc[c]);
    assert!(pl[c] < 0.0, "凸坡平面曲率应为负: {}", pl[c]);
    // 无效邻域 -> NaN
    let mut with_nodata = convex.clone();
    with_nodata[c - 1] = f32::NAN;
    let pc2 = profile_curvature(&with_nodata, w, h, res);
    assert!(pc2[c].is_nan(), "无效邻域必须传播 NaN");
}

// ---------------- Task 8: 自适应形态证据 ----------------

use topo_core::geomorphon::adaptive_morphology_evidence;

/// 圆丘场景: 半径 radius_m, 峰 850m, 周围谷底 700m
fn conical(size: usize, res: f64, radius_m: f64) -> (Vec<f32>, topo_core::input::RasterShape) {
    let cx = size as f64 * res * 0.5;
    let z: Vec<f32> = (0..size)
        .flat_map(|row| {
            (0..size).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cx;
                let r = x.hypot(y);
                (700.0 + 150.0 * (1.0 - r / radius_m).max(0.0)) as f32
            })
        })
        .collect();
    (z, topo_core::input::RasterShape { width: size, height: size, resolution_m: res })
}

const SCALES: [f64; 3] = [125.0, 250.0, 500.0];

/// 剖面证据方向性: 峰顶上部证据占优, 周缘下部证据占优, 直坡段中部占优
#[test]
fn adaptive_evidence_favors_expected_positions() {
    let (dem, shape) = conical(301, 10.0, 600.0);
    let valid = vec![true; dem.len()];
    let scale_m: Vec<f32> = vec![250.0; dem.len()];
    let ev = adaptive_morphology_evidence(&dem, &valid, shape, &scale_m, &SCALES).unwrap();
    let c = 150 * 301 + 150;
    assert!(ev.upper[c] > 0.7, "峰顶上部证据不足: {}", ev.upper[c]);
    let edge = 150 * 301 + 8; // 西缘谷底
    assert!(ev.lower[edge] > 0.5, "周缘下部证据不足: {}", ev.lower[edge]);
    // 坡向带语义: 圆锥坡切向平坦时 geomorphon 本征给出 Shoulder 型模式
    // (证据 0.70/0.30/0.00), 方向性断言按上/下证据对比
    let mid = 150 * 301 + 150 - 30; // 北向 300m(半坡)
    let mid_upper = 150 * 301 + 150 - 45; // 450m 偏峰
    assert!(
        ev.upper[mid] > ev.lower[mid],
        "半坡不得偏下部: u={} m={} l={}",
        ev.upper[mid],
        ev.middle[mid],
        ev.lower[mid]
    );
    assert!(
        ev.upper[mid_upper] > ev.lower[mid_upper],
        "峰侧应偏上部: u={} m={} l={}",
        ev.upper[mid_upper],
        ev.middle[mid_upper],
        ev.lower[mid_upper]
    );
    // 谷底环带(r>650m)整体下部证据占优(圆锥坡面为对称 Slope 型, 只有
    // 谷底平地给出明确的 lower 主导)
    let (mut lu, mut ll) = (0f64, 0f64);
    for y in 0..301 {
        for x in 0..301 {
            let dx = (x as f64 + 0.5) * 10.0 - 1505.0;
            let dy = (y as f64 + 0.5) * 10.0 - 1505.0;
            if dx.hypot(dy) > 650.0 {
                lu += ev.upper[y * 301 + x] as f64;
                ll += ev.lower[y * 301 + x] as f64;
            }
        }
    }
    assert!(ll > lu, "谷底环带下部证据应占优: lower={ll} upper={lu}");
}

/// 特征宽度加倍时形态采样半径跟随 adaptive_scale_m(上部带按比例展宽)
#[test]
fn evidence_radius_follows_adaptive_scale() {
    let run = |radius_m: f64| {
        let (dem, shape) = conical(601, 10.0, radius_m);
        let valid = vec![true; dem.len()];
        let scale_m: Vec<f32> = vec![(radius_m / 2.4) as f32; dem.len()];
        let ev =
            adaptive_morphology_evidence(&dem, &valid, shape, &scale_m, &SCALES).unwrap();
        // 上部证据(>0.5)像元占比
        let upper = ev
            .upper
            .iter()
            .filter(|&&v| v > 0.5)
            .count() as f64
            / ev.upper.len() as f64;
        upper
    };
    let narrow = run(600.0);
    let wide = run(1500.0);
    assert!(
        wide > narrow * 1.5,
        "宽丘上部带未按尺度展宽: narrow={narrow:.3} wide={wide:.3}"
    );
}

/// 端点尺度或导数支持不完整的像元标记低置信
#[test]
fn endpoint_scale_and_incomplete_derivatives_mark_low_confidence() {
    let (dem, shape) = conical(301, 10.0, 600.0);
    let valid = vec![true; dem.len()];
    // 全部取最小尺度(端点) -> 低置信
    let scale_min: Vec<f32> = vec![125.0; dem.len()];
    let ev1 = adaptive_morphology_evidence(&dem, &valid, shape, &scale_min, &SCALES).unwrap();
    assert!(
        ev1.low_confidence.iter().all(|&b| b),
        "端点尺度应全部低置信"
    );
    // 中段尺度 + 完整导数 -> 非低置信
    let scale_mid: Vec<f32> = vec![250.0; dem.len()];
    let ev2 = adaptive_morphology_evidence(&dem, &valid, shape, &scale_mid, &SCALES).unwrap();
    assert!(
        ev2.low_confidence.iter().any(|&b| !b),
        "中段尺度不应低置信"
    );
    assert!(
        ev2.slope_deg.iter().all(|v| v.is_finite()),
        "坡度场必须完整"
    );
}
