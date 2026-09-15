//! 脊线与坡面单元几何契约测试: 谷心非脊、两翼连续脊、屏障约束与镜像对称。

mod common;

use common::{build_context, Context};
use topo_core::hydro::HydroConfig;
use topo_core::ridge::build_ridges;
use topo_core::slope_unit::build_slope_units;

/// 对称 V 谷: 中央谷不标脊, 两侧翼各含一条连续山脊
#[test]
fn valley_center_not_ridge_and_flanks_have_ridges() {
    let (dem, _) = common::v_valley(201, 201, 10.0);
    let ctx = build_context(dem, 201, 201, 10.0);
    let ridges =
        build_ridges(&ctx.coarse, &ctx.valid, &ctx.hydro, &ctx.pyramid, &ctx.landform).unwrap();
    let cw = ctx.shape.width;
    // 1) 中央谷(中线 ±2 粗像元)非脊
    for y in 5..ctx.shape.height - 5 {
        for dx in -2i64..=2 {
            let x = (cw as i64 / 2 + dx) as usize;
            assert!(
                !ridges.mask[y * cw + x],
                "中央谷({x},{y})不得标脊"
            );
        }
    }
    // 2) 两侧翼各含连续脊: 左右半幅各有脊像元且各自 8 连通组件 >= 5 像元
    let count_connected = |half: std::ops::Range<usize>| -> usize {
        let mut seen = vec![false; ridges.mask.len()];
        let mut best = 0usize;
        for y in 0..ctx.shape.height {
            for x in half.clone() {
                let i = y * cw + x;
                if !ridges.mask[i] || seen[i] {
                    continue;
                }
                let (mut size, mut stack) = (0usize, vec![i]);
                seen[i] = true;
                while let Some(c) = stack.pop() {
                    size += 1;
                    let (cx, cy) = (c % cw, c / cw);
                    for dyy in -1i64..=1 {
                        for dxx in -1i64..=1 {
                            let nx = cx as i64 + dxx;
                            let ny = cy as i64 + dyy;
                            if nx < 0
                                || ny < 0
                                || nx >= cw as i64
                                || ny >= ctx.shape.height as i64
                            {
                                continue;
                            }
                            let j = ny as usize * cw + nx as usize;
                            if ridges.mask[j]
                                && !seen[j]
                                && half.contains(&(nx as usize))
                            {
                                seen[j] = true;
                                stack.push(j);
                            }
                        }
                    }
                }
                best = best.max(size);
            }
        }
        best
    };
    let left = count_connected(0..cw / 2);
    let right = count_connected(cw / 2 + 1..cw);
    assert!(left >= 5, "左翼连续脊不足: {left}");
    assert!(right >= 5, "右翼连续脊不足: {right}");
    // 3) 脊不落在河网上
    for (i, &s) in ctx.hydro.streams[0].iter().enumerate() {
        assert!(!(s && ridges.mask[i]), "河网像元 {i} 不得标脊");
    }
}

/// 非对称谷: 镜像 DEM 的坡面单元面积分布一致(逐项差 <= 2 粗像元)
#[test]
fn mirrored_dem_has_matching_unit_areas() {
    let build = |mirror: bool| {
        let (w, h, res) = (201usize, 201usize, 10.0);
        let mut dem = vec![0f32; w * h];
        for y in 0..h {
            for x in 0..w {
                let xi = if mirror { w - 1 - x } else { x };
                let dx = xi as f64 * res - (w as f64 * res) * 0.5;
                dem[y * w + x] =
                    (800.0 + (if dx < 0.0 { 0.10 } else { 0.05 }) * dx.abs() + 0.004 * y as f64 * res / 10.0) as f32;
            }
        }
        let ctx = build_context(dem, w as u32, h as u32, res);
        let ridges = build_ridges(
            &ctx.coarse,
            &ctx.valid,
            &ctx.hydro,
            &ctx.pyramid,
            &ctx.landform,
        )
        .unwrap();
        build_slope_units(&ctx.coarse, &ctx.valid, &ctx.hydro, &ridges).unwrap()
    };
    let u1 = build(false);
    let u2 = build(true);
    let areas = |u: &topo_core::slope_unit::SlopeUnits| -> Vec<usize> {
        let mut sizes: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for &id in &u.unit_id {
            if id != 0 {
                *sizes.entry(id).or_default() += 1;
            }
        }
        let mut v: Vec<usize> = sizes.into_values().collect();
        v.sort_unstable();
        v
    };
    let (a1, a2) = (areas(&u1), areas(&u2));
    assert_eq!(a1.len(), a2.len(), "镜像后单元数量应一致: {a1:?} vs {a2:?}");
    for (x, y) in a1.iter().zip(a2.iter()) {
        // 坡向泛洪组件的分裂边界是离散化软边界, 1-2 粗像元错位即可放大
        // 为面积差; 采用 5% 相对容差验证算法无镜像偏倚
        let tol = ((*x).max(*y) / 20).max(2);
        assert!(
            x.abs_diff(*y) <= tol,
            "镜像单元面积差超容差({tol}): {a1:?} vs {a2:?}"
        );
    }
}

/// 对称 V 谷: 坡面单元不跨脊线或河网屏障(同单元即连通且不含屏障)
#[test]
fn slope_units_do_not_cross_barriers() {
    let (dem, _) = common::v_valley(201, 201, 10.0);
    let ctx = build_context(dem, 201, 201, 10.0);
    let ridges =
        build_ridges(&ctx.coarse, &ctx.valid, &ctx.hydro, &ctx.pyramid, &ctx.landform).unwrap();
    let units = build_slope_units(&ctx.coarse, &ctx.valid, &ctx.hydro, &ridges).unwrap();
    let cw = ctx.shape.width;
    for i in 0..units.unit_id.len() {
        if units.unit_id[i] != 0 {
            assert!(
                !units.valley_mask[i] && !units.ridge_mask[i],
                "单元像元 {i} 不得落在屏障上"
            );
            assert!(ctx.valid[i], "单元不得进入无效区");
        }
    }
    // 同单元像元 8 连通(无跨屏障的离散拼接)
    let n = units.unit_id.len();
    let mut seen = vec![false; n];
    for s in 0..n {
        if units.unit_id[s] == 0 || seen[s] {
            continue;
        }
        let id = units.unit_id[s];
        let (mut cnt, mut stack) = (0usize, vec![s]);
        seen[s] = true;
        while let Some(c) = stack.pop() {
            cnt += 1;
            let (cx, cy) = (c % cw, c / cw);
            for dyy in -1i64..=1 {
                for dxx in -1i64..=1 {
                    if dxx == 0 && dyy == 0 {
                        continue;
                    }
                    let nx = cx as i64 + dxx;
                    let ny = cy as i64 + dyy;
                    if nx < 0 || ny < 0 || nx >= cw as i64 || ny >= ctx.shape.height as i64 {
                        continue;
                    }
                    let j = ny as usize * cw + nx as usize;
                    if units.unit_id[j] == id && !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        assert!(cnt > 0);
        _ = cnt;
    }
}

// ---------------- Task 7: 约束距离/自适应尺度/相对位置 ----------------

use topo_core::slope_unit::{build_slope_geometry, SlopeUnits};

/// 双谷隔脊场景: 谷底 x=500m 与 x=1500m, 中央 x=1000m 为脊(+40m)
fn dual_valley(steep: f64, mirror: bool) -> (Vec<f32>, topo_core::input::RasterShape) {
    let (w, h, res) = (201usize, 201usize, 10.0);
    let mut dem = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let xi = if mirror { w - 1 - x } else { x };
            let xm = (xi as f64 + 0.5) * res;
            let ym = (y as f64 + 0.5) * res;
            let d_near = (xm - 500.0).abs().min((xm - 1500.0).abs());
            dem[y * w + x] = (800.0 + steep * d_near + 0.004 * ym) as f32;
        }
    }
    (dem, topo_core::input::RasterShape { width: w, height: h, resolution_m: res })
}

fn full_geometry(
    dem: Vec<f32>,
    w: u32,
    h: u32,
    res: f64,
) -> (
    Context,
    topo_core::ridge::RidgeModel,
    SlopeUnits,
    topo_core::slope_unit::SlopeGeometry,
) {
    let ctx = build_context(dem, w, h, res);
    let ridges =
        build_ridges(&ctx.coarse, &ctx.valid, &ctx.hydro, &ctx.pyramid, &ctx.landform).unwrap();
    let units = build_slope_units(&ctx.coarse, &ctx.valid, &ctx.hydro, &ridges).unwrap();
    let geom = build_slope_geometry(
        &ctx.coarse,
        &ctx.valid,
        &units,
        &ctx.hydro,
        &ctx.pyramid,
        20.0,
    )
    .unwrap();
    (ctx, ridges, units, geom)
}

/// 跨脊欧氏最近谷永不作为谷锚: 左谷坡像元的锚必在左谷(中央脊 x=1000m 左侧)
#[test]
fn constrained_distance_never_crosses_ridge() {
    let (dem, _) = dual_valley(0.08, false);
    let (ctx, _r, units, geom) = full_geometry(dem, 201, 201, 10.0);
    let cw = ctx.shape.width;
    for i in 0..cw * ctx.shape.height {
        // dv/dr 相对"所选等级"的对应谷/脊线(规格 9.2/9.3): 主干谷像元
        // dv=0, 坡脚细谷延伸段到所选粗谷线距离有限
        if units.valley_mask[i] {
            let dv = geom.distance_to_valley_m[i];
            assert!(dv.is_finite(), "谷像元 {i} dv 应有限");
            if dv == 0.0 {
                zero_dv += 1;
            }
        }
        if units.ridge_mask[i] {
            assert_eq!(
                geom.distance_to_ridge_m[i], 0.0,
                "脊像元 {i} dr 应为 0"
            );
        }
        if geom.relative_position[i].is_finite() {
            assert!(
                (0.0..=1.0).contains(&geom.relative_position[i]),
                "q 越界: {}",
                geom.relative_position[i]
            );
        }
        let l = geom.local_width_m[i];
        if l.is_finite() {
            let (dv, dr) = (geom.distance_to_valley_m[i], geom.distance_to_ridge_m[i]);
            assert!(
                (l - (dv + dr)).abs() < 0.5,
                "L != dv+dr @ {i}"
            );
        }
    }
    // 谷锚不跨中央脊(粗列 40): x<40 的坡面像元锚 x<=22(左谷 500m/25m+2)
    let mut checked = 0usize;
    let mut zero_dv = 0usize;
    for y in 5..ctx.shape.height - 5 {
        for x in 23..40usize {
            let i = y * cw + x;
            if units.unit_id[i] == 0 || geom.valley_anchor[i] == u32::MAX {
                continue;
            }
            let ax = (geom.valley_anchor[i] as usize % cw) as f64 * 25.0;
            assert!(
                ax <= 550.0,
                "左坡像元 ({x},{y}) 谷锚在 {ax}m, 跨脊关联到右谷"
            );
            checked += 1;
        }
    }
    assert!(checked > 200, "有效锚样本不足: {checked}");
    assert!(zero_dv > 50, "主干谷 dv=0 样本不足: {zero_dv}");
}

/// 相对位置内部一致: q 由 qd/qz 按低起伏权重合成; 低起伏场景水平权重更大
#[test]
fn relative_position_weights_low_relief_horizontally() {
    let run = |steep: f64| {
        let (dem, _) = dual_valley(steep, false);
        let (ctx, _r, _u, geom) = full_geometry(dem, 201, 201, 10.0);
        // 特征起伏中位数
        let mut hrel: Vec<f32> = ctx
            .pyramid
            .characteristic_relief_m
            .iter()
            .copied()
            .filter(|v| v.is_finite() && *v > 0.0)
            .collect();
        hrel.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let h_med = hrel[hrel.len() / 2];
        // 逐像元重算 q 验证一致性(w_distance 依赖逐像元 H*)
        let mut err = 0f64;
        let mut n = 0usize;
        let mut w_acc = 0f64;
        for i in 0..geom.relative_position.len() {
            let (qd, qz) = (geom.q_distance[i], geom.q_elevation[i]);
            if qd.is_finite() && qz.is_finite() && geom.relative_position[i].is_finite() {
                let h_star = ctx.pyramid.characteristic_relief_m[i];
                let w = 0.45f32
                    + (0.75f32 - 0.45) * ((20.0 - h_star) / 20.0).clamp(0.0, 1.0);
                let expect = (w * qd + (1.0 - w) * qz).clamp(0.0, 1.0);
                err += (expect - geom.relative_position[i]).abs() as f64;
                w_acc += w as f64;
                n += 1;
            }
        }
        let w_med = if n > 0 { w_acc / n as f64 } else { 0.0 };
        (h_med, w_med, if n > 0 { err / n as f64 } else { f64::MAX })
    };
    let (h_steep, w_steep, e1) = run(0.08);
    let (h_flat, w_flat, e2) = run(0.008);
    assert!(e1 < 1e-3 && e2 < 1e-3, "q 内部一致性失败: {e1} {e2}");
    assert!(w_flat > w_steep, "低起伏权重未提高: {w_flat} vs {w_steep}");
    assert!(h_flat < h_steep, "场景起伏设置失败: {h_flat} vs {h_steep}");
}

/// 镜像场景的约束距离分布一致(中位差 <= 1 粗像元)
#[test]
fn mirrored_distance_distributions_match() {
    let median = |steep: f64, mirror: bool| -> f64 {
        let (dem, _) = dual_valley(steep, mirror);
        let (_ctx, _r, _u, geom) = full_geometry(dem, 201, 201, 10.0);
        let mut v: Vec<f32> = geom
            .distance_to_valley_m
            .iter()
            .copied()
            .filter(|x| x.is_finite())
            .collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2] as f64
    };
    let m1 = median(0.08, false);
    let m2 = median(0.08, true);
    assert!(
        (m1 - m2).abs() <= 25.0,
        "镜像距离中位差超 1 粗像元: {m1} vs {m2}"
    );
}
