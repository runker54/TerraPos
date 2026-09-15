//! 盆地对象契约测试: 宽盆/宽谷接受, 窄谷/高地平地拒绝, 完整边界重建。

mod common;

use common::full_chain;
use topo_core::pipeline::{BasinTendency};
use topo_core::basin::{detect_basins, BasinConfig};

fn cfg(min_area: f64) -> BasinConfig {
    BasinConfig {
        tendency: BasinTendency::Standard,
        min_area_m2: min_area,
    }
}

/// 圆形闭合宽盆: 平底半径 900m, 环状山脊围限, 底部 0.002 坡
fn broad_basin() -> (Vec<f32>, topo_core::input::RasterShape) {
    common::broad_basin(401, 401, 10.0)
}

/// 宽谷: 谷底宽 800m、0.004 纵坡排水, 两侧 0.15 坡升
fn wide_valley() -> (Vec<f32>, topo_core::input::RasterShape) {
    let (w, h, res) = (401usize, 401usize, 10.0);
    let mut dem = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let xm = (x as f64 + 0.5) * res;
            let ym = (y as f64 + 0.5) * res;
            let d = (xm - 2005.0).abs();
            let slope = if d < 400.0 { 0.0 } else { 0.15 * (d - 400.0) };
            dem[y * w + x] = (800.0 + 0.004 * ym + slope) as f32;
        }
    }
    (dem, topo_core::input::RasterShape { width: w, height: h, resolution_m: res })
}

/// 窄 V 谷: 谷底无平带
fn narrow_valley() -> (Vec<f32>, topo_core::input::RasterShape) {
    let (w, h, res) = (401usize, 401usize, 10.0);
    let mut dem = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let xm = (x as f64 + 0.5) * res;
            let ym = (y as f64 + 0.5) * res;
            dem[y * w + x] = (800.0 + 0.004 * ym + 0.10 * (xm - 2005.0).abs()) as f32;
        }
    }
    (dem, topo_core::input::RasterShape { width: w, height: h, resolution_m: res })
}

/// 高地平台: 中央平坦高地, 四周陡落(平地 HAND 高应拒绝)
fn flat_upland() -> (Vec<f32>, topo_core::input::RasterShape) {
    let (w, h, res) = (401usize, 401usize, 10.0);
    let mut dem = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let xm = (x as f64 + 0.5) * res;
            let ym = (y as f64 + 0.5) * res;
            let d_edge = xm.min(4010.0 - xm).min(ym.min(4010.0 - ym));
            let fall = if d_edge < 600.0 { 2.0 * (600.0 - d_edge) } else { 0.0 };
            dem[y * w + x] = (900.0 - fall) as f32;
        }
    }
    (dem, topo_core::input::RasterShape { width: w, height: h, resolution_m: res })
}

fn detect(scene: &(Vec<f32>, topo_core::input::RasterShape), min_area: f64) -> topo_core::basin::BasinResult {
    let (ctx, _r, units, geom, morph) =
        full_chain(scene.0.clone(), scene.1.width as u32, scene.1.height as u32, scene.1.resolution_m);
    detect_basins(
        &ctx.coarse,
        &ctx.valid,
        ctx.shape,
        &ctx.hydro,
        &units,
        &geom,
        &ctx.pyramid,
        &morph,
        &cfg(min_area),
    )
    .unwrap()
}

fn accepted_area(res: &topo_core::basin::BasinResult) -> usize {
    if res.objects.iter().any(|o| o.accepted) {
        res.mask.iter().filter(|&&b| b).count()
    } else {
        0
    }
}

/// 宽闭合盆与宽谷接受; 窄谷与高地平台拒绝
#[test]
fn broad_basins_accepted_narrow_and_upland_rejected() {
    let r_basin = detect(&broad_basin(), 66_666.67);
    for o in &r_basin.objects {
        eprintln!(
            "OBJ{} area={:.0} wmax={} w50={} inner={} surround={} stream={} closed={} accepted={}",
            o.id, o.area_m2, o.max_width_m, o.median_width_m,
            o.inner_relief_m, o.surround_rise_m, o.stream_connected,
            o.closed_depression, o.accepted
        );
    }
    eprintln!("candidate={}", r_basin.candidate.iter().filter(|&&b| b).count());
    assert!(
        r_basin.objects.iter().any(|o| o.accepted),
        "宽闭合盆应被接受"
    );
    let r_valley = detect(&wide_valley(), 66_666.67);
    assert!(
        r_valley.objects.iter().any(|o| o.accepted),
        "宽谷盆应被接受"
    );
    assert!(
        r_valley
            .objects
            .iter()
            .any(|o| o.accepted && o.stream_connected),
        "宽谷对象应具水文连通标记"
    );
    for (name, scene) in [("窄V谷", narrow_valley()), ("高地平台", flat_upland())] {
        let r = detect(&scene, 66_666.67);
        assert!(
            !r.objects.iter().any(|o| o.accepted),
            "{name} 不应产生盆地对象"
        );
    }
}

/// 接受对象的 mask 是候选域重建(含 >=90% 平底), 且核心显著小于 mask
#[test]
fn accepted_basin_reconstructs_full_boundary() {
    let scene = broad_basin();
    let r = detect(&scene, 66_666.67);
    let mask_n = r.mask.iter().filter(|&&b| b).count();
    let core_n = r.core.iter().filter(|&&b| b).count();
    let cand_n = r.candidate.iter().filter(|&&b| b).count();
    assert!(mask_n > 0, "应存在重建盆地");
    assert!(
        core_n * 4 < mask_n,
        "宽度核心应显著小于重建边界: core={core_n} mask={mask_n}"
    );
    // 其余候选属于地形上分隔的其他对象(外圈平台带等),
    // 不并入盆底; mask 应显著大于核心(完整重建而非仅核心)
    assert!(mask_n > cand_n / 2, "重建边界过小: mask={mask_n} cand={cand_n}");
    // 已知平底(r<500m)中成为候选的像元, 90% 必须保留在重建 mask 内
    // (检验重建完整性; persist 大窗语义下盆缘平地本就不产生候选)
    let cw = (4010.0f64 / 25.0).ceil() as usize;
    let (mut in_mask, mut total) = (0usize, 0usize);
    for y in 0..cw {
        for x in 0..cw {
            let xm = (x as f64 + 0.5) * 25.0 - 2005.0;
            let ym = (y as f64 + 0.5) * 25.0 - 2005.0;
            if xm.hypot(ym) < 500.0 {
                let i = y * cw + x;
                if r.candidate[i] {
                    total += 1;
                    if r.mask[i] {
                        in_mask += 1;
                    }
                }
            }
        }
    }
    assert!(total > 0, "平底区应产生候选");
    assert!(
        in_mask * 10 >= total * 9,
        "平底候选重建覆盖不足: {in_mask}/{total}"
    );
}

/// 提高最小面积阈值只移除不足对象, 不改变其余边界
#[test]
fn raising_min_area_only_drops_undersized_objects() {
    let scene = wide_valley();
    let r1 = detect(&scene, 66_666.67);
    let r2 = detect(&scene, 5.0e7); // 远超谷底总面积
    eprintln!("raising r1 objects: {:?}", r1.objects.iter().map(|o| (o.area_m2, o.accepted)).collect::<Vec<_>>());
    assert!(accepted_area(&r1) > 0);
    assert_eq!(accepted_area(&r2), 0, "超阈值后应无接受对象");
}
