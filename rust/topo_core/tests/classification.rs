//! 分类契约测试
//!
//! 使用共享合成夹具。: 模糊隶属度、置信度与确定性坡位(Task 11 后扩展组合与清理)。

mod common;

// ---------------- Task 9: 模糊上中下隶属度 ----------------

use topo_core::geomorphon::{Landform, MorphEvidence};
use topo_core::slope_position::{classify_slope_positions, SlopePosition};

fn neutral_morph(n: usize) -> MorphEvidence {
    MorphEvidence {
        form: vec![Landform::Flat; n],
        upper: vec![0.0; n],
        middle: vec![0.0; n],
        lower: vec![0.0; n],
        slope_deg: vec![5.0; n],
        profile_curvature: vec![0.0; n],
        plan_curvature: vec![0.0; n],
        low_confidence: vec![false; n],
    }
}

/// q=0.1/0.5/0.9 处隶属度排序正确; 中心平局取中部
#[test]
fn fuzzy_membership_orders_by_relative_position() {
    let qs: Vec<f32> = vec![0.1, 0.5, 0.9];
    let n = qs.len();
    let morph = neutral_morph(n);
    let m = classify_slope_positions(&qs, &morph, &vec![true; n]).unwrap();
    assert_eq!(m.raw[0], SlopePosition::Lower, "q=0.1 应为下部");
    assert_eq!(m.raw[1], SlopePosition::Middle, "q=0.5 中心应为中部(平局取中)");
    assert_eq!(m.raw[2], SlopePosition::Upper, "q=0.9 应为上部");
}

/// 0.38/0.62 过渡带两侧隶属度连续(无跳变)
#[test]
fn fuzzy_membership_is_continuous_at_transitions() {
    let qs: Vec<f32> = vec![0.37, 0.39, 0.61, 0.63];
    let morph = neutral_morph(qs.len());
    let m = classify_slope_positions(&qs, &morph, &vec![true; qs.len()]).unwrap();
    for k in [0usize, 2] {
        let d_upper = (m.upper[k] - m.upper[k + 1]).abs();
        let d_lower = (m.lower[k] - m.lower[k + 1]).abs();
        assert!(d_upper < 0.15, "上隶属度跨 0.62 跳变: {d_upper}");
        assert!(d_lower < 0.15, "下隶属度跨 0.38 跳变: {d_lower}");
    }
}

/// 形态证据可扭转接近的平局: 中部位置 + 强峰证据 -> 上部
#[test]
fn morphology_resolves_near_tie() {
    let n = 1;
    let mut morph = neutral_morph(n);
    morph.upper[0] = 5.0; // 强峰证据
    let m = classify_slope_positions(&[0.5f32], &morph, &[true]).unwrap();
    assert_eq!(m.raw[0], SlopePosition::Upper, "强峰证据应扭转中部主导");
}

/// 置信度 = 最大隶属度 - 次大隶属度; 低置信折扣 0.6 不改变排序
#[test]
fn confidence_and_low_confidence_discount() {
    let n = 2;
    let qs = [0.9f32, 0.9];
    let mut morph = neutral_morph(n);
    let m = classify_slope_positions(&qs, &morph, &vec![true; n]).unwrap();
    let mut s: Vec<f32> = vec![m.upper[0], m.middle[0], m.lower[0]];
    s.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let expect = s[0] - s[1];
    assert!((m.confidence[0] - expect).abs() < 1e-4, "置信度错误");
    // 低置信(端点尺度) -> 乘 0.6
    morph.low_confidence[1] = true;
    let m2 = classify_slope_positions(&qs, &morph, &vec![true; n]).unwrap();
    assert!((m2.confidence[1] - m.confidence[1] * 0.6).abs() < 1e-4);
    assert_eq!(m2.raw[1], m.raw[0], "折扣不得改变排序");
}

/// q 无效(锚点缺失)时仅凭形态证据分类且置信度打折
#[test]
fn missing_geometry_falls_back_to_morphology() {
    let n = 1;
    let mut morph = neutral_morph(n);
    morph.lower[0] = 1.0;
    morph.form[0] = Landform::Pit;
    let m = classify_slope_positions(&[f32::NAN], &morph, &[true]).unwrap();
    assert_eq!(m.raw[0], SlopePosition::Lower, "仅形态证据时应为下部");
    assert!(m.confidence[0] <= 0.6, "几何缺失置信度应打折");
}

// ---------------- Task 11: 编码组合与受约束清理 ----------------

use common::full_chain;
use topo_core::basin::{detect_basins, BasinConfig};
use topo_core::pipeline::{BasinTendency};
use topo_core::postprocess::{compose_codes, constrained_cleanup};

/// 组合规则: 盆地优先=1; 丘陵(<500m) 下/中/上=5/4/3; 山地 8/7/6;
/// 无效=0; 保留码 2 永不出现
#[test]
fn compose_codes_basin_overrides_and_elevation_split() {
    let w = 8usize;
    let h = 4usize;
    let n = w * h;
    let mut elevation = vec![400.0f32; n]; // 上半幅丘陵
    for y in 2..4 {
        for x in 0..w {
            elevation[y * w + x] = 600.0; // 下半幅山地
        }
    }
    let valid = vec![true; n];
    let basin = topo_core::basin::BasinResult {
        candidate: vec![false; n],
        core: vec![false; n],
        mask: {
            let mut m = vec![false; n];
            m[0] = true;
            m[1] = true;
            m
        },
        objects: vec![],
    };
    let mut positions = topo_core::slope_position::Memberships {
        upper: vec![0.0; n],
        middle: vec![0.0; n],
        lower: vec![0.0; n],
        raw: vec![SlopePosition::Lower; n],
        confidence: vec![0.5; n],
    };
    // 上: index 2 丘陵 Upper / index 2+w 山地 Lower / index 2+2w 山地 Middle
    positions.raw[2] = SlopePosition::Upper;
    positions.raw[2 + 2 * w] = SlopePosition::Lower;
    positions.raw[2 + 3 * w] = SlopePosition::Middle;
    let fc = compose_codes(&elevation, &valid, &basin, &positions, 500.0).unwrap();
    // 盆地覆盖(index 0/1)
    assert_eq!(fc.terrain[0], 1);
    assert_eq!(fc.terrain[1], 1);
    // 丘陵(400m): 上3
    assert_eq!(fc.terrain[2], 3);
    // 山地(600m): 下8 / 中7
    assert_eq!(fc.terrain[2 + 2 * w], 8);
    assert_eq!(fc.terrain[2 + 3 * w], 7);
    // 保留码
    assert!(!fc.terrain.contains(&2), "保留码 2 不得出现");
    // 亚类独立: 山地行亚类=低山(3), 丘陵行=低丘(1)
    assert_eq!(fc.geomorph_subclass[2 + 2 * w], 3);
    assert_eq!(fc.geomorph_subclass[2], 1);
    // NoData
    let mut valid2 = valid.clone();
    valid2[5] = false;
    let fc2 = compose_codes(&elevation, &valid2, &basin, &positions, 500.0).unwrap();
    assert_eq!(fc2.terrain[5], 0);
}

/// 单像元斑贴着山脊时, 清理不得把它搬到脊对侧
#[test]
fn cleanup_never_crosses_ridge_barrier() {
    let (dem, _) = common::v_valley(201, 201, 10.0);
    let (ctx, ridges, units, geom, morph) = full_chain(dem, 201, 201, 10.0);
    let basin = detect_basins(
        &ctx.coarse,
        &ctx.valid,
        ctx.shape,
        &ctx.hydro,
        &units,
        &geom,
        &ctx.pyramid,
        &morph,
        &BasinConfig { tendency: BasinTendency::Standard, min_area_m2: 66_666.67 },
    )
    .unwrap();
    let mut fc = topo_core::postprocess::FinalClassification {
        terrain: vec![6; ctx.coarse.len()],
        geomorph_subclass: vec![3; ctx.coarse.len()],
        confidence: vec![0.9; ctx.coarse.len()],
    };
    // 在中央谷(两侧山脊之间)放一个 5 类单像元斑
    let cw = ctx.shape.width;
    let mid = cw / 2;
    let spot = 20 * cw + mid;
    fc.terrain[spot] = 5;
    let before_ridge = fc.terrain.clone();
    constrained_cleanup(&mut fc, &morph_low(&morph, ctx.coarse.len()), &units, &geom, &basin, ctx.shape, 1.0).unwrap();
    // 斑要么保留要么并到同侧相邻类, 山脊另一侧像元不得被改成 5
    for y in 0..ctx.shape.height {
        for x in 0..cw {
            let i = y * cw + x;
            if before_ridge[i] == 6 {
                let cross_ridge = (0..cw).any(|x2| {
                    let between = ridges.mask[(y * cw + ((x + x2) / 2))] && units.ridge_mask[y * cw + ((x + x2) / 2)];
                    between && fc.terrain[i] == 5
                });
                assert!(!cross_ridge, "像元 {i} 跨脊被清理成 5");
            }
        }
    }
}

fn morph_low(
    morph: &topo_core::geomorphon::MorphEvidence,
    n: usize,
) -> topo_core::slope_position::Memberships {
    let q: Vec<f32> = vec![0.5; n];
    let m = topo_core::slope_position::classify_slope_positions(&q, morph, &vec![true; n]).unwrap();
    m
}

/// 500 条确定性脊→谷轨迹在受约束清理后 >=95% 满足单调顺序
#[test]
fn traces_monotonic_after_constrained_cleanup() {
    let (dem, _) = common::nested_ridges(401, 401, 10.0);
    let (ctx, ridges, units, geom, morph) = full_chain(dem, 401, 401, 10.0);
    let basin = detect_basins(
        &ctx.coarse,
        &ctx.valid,
        ctx.shape,
        &ctx.hydro,
        &units,
        &geom,
        &ctx.pyramid,
        &morph,
        &BasinConfig { tendency: BasinTendency::Standard, min_area_m2: 66_666.67 },
    )
    .unwrap();
    let positions =
        topo_core::slope_position::classify_slope_positions(&geom.relative_position, &morph, &ctx.valid)
            .unwrap();
    let mut fc = compose_codes(&ctx.coarse, &ctx.valid, &basin, &positions, 500.0).unwrap();
    constrained_cleanup(&mut fc, &positions, &units, &geom, &basin, ctx.shape, 1.0).unwrap();
    let cw = ctx.shape.width;
    let rank = |c: u8| match c {
        3 | 6 => 2i32,
        4 | 7 => 1,
        5 | 8 | 1 => 0,
        _ => -1,
    };
    // 轨迹: 从每个脊像元沿 dv 递减步进至谷锚/无进展
    let mut traces = Vec::new();
    'seeds: for s in 0..ctx.coarse.len() {
        if !ridges.mask[s] || !ctx.valid[s] {
            continue;
        }
        if traces.len() >= 500 {
            break 'seeds;
        }
        let mut trace = vec![s];
        let mut cur = s;
        loop {
            let x = (cur % cw) as i64;
            let y = (cur / cw) as i64;
            let mut next = None;
            let mut best_dv = geom.distance_to_valley_m[cur];
            for (dx, dy) in [
                (-1i64, 0i64),
                (1, 0),
                (0, -1),
                (0, 1),
                (-1, -1),
                (1, -1),
                (-1, 1),
                (1, 1),
            ] {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= cw as i64 || ny >= ctx.shape.height as i64 {
                    continue;
                }
                let j = (ny as usize) * cw + nx as usize;
                if !ctx.valid[j] || trace.contains(&j) {
                    continue;
                }
                let dv = geom.distance_to_valley_m[j];
                if dv.is_finite() && dv < best_dv {
                    best_dv = dv;
                    next = Some(j);
                }
            }
            match next {
                Some(j) => {
                    trace.push(j);
                    cur = j;
                }
                None => break,
            }
        }
        if trace.len() >= 5 {
            traces.push(trace);
        }
    }
    assert!(traces.len() >= 500, "轨迹不足 500: {}", traces.len());
    let mut monotonic = 0usize;
    for trace in &traces {
        let mut prev = 2i32;
        let mut ok = true;
        for &i in trace {
            let r = rank(fc.terrain[i]);
            if r < 0 {
                continue;
            }
            if r > prev {
                ok = false;
                break;
            }
            prev = r;
        }
        if ok {
            monotonic += 1;
        }
    }
    let rate = monotonic as f64 / traces.len() as f64;
    assert!(rate >= 0.95, "单调率 {rate:.3} < 95%");
}
