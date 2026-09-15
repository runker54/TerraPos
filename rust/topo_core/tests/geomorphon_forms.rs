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
