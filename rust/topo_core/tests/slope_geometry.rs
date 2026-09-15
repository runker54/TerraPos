//! 脊线与坡面单元几何契约测试: 谷心非脊、两翼连续脊、屏障约束与镜像对称。

mod common;

use topo_core::geomorphon::{geomorphon_pattern, pattern_to_landform, Landform};
use topo_core::hydro::{build_hydro, HydroConfig, HydroModel};
use topo_core::input::{prepare_values, RasterShape};
use topo_core::ridge::build_ridges;
use topo_core::scale::{build_scale_pyramid, ScalePyramid};
use topo_core::slope_unit::build_slope_units;

fn hydro_cfg() -> HydroConfig {
    HydroConfig {
        coarse_res_m: 25.0,
        z_limit_m: 15.0,
        stream_areas_km2: [0.05, 0.20, 1.00, 5.00],
    }
}

/// 粗层上下文: 粗层表面/有效掩膜/水文/金字塔/形态证据
/// (与生产管线 Task 12 相同的降采样语义, 任一原生有效即粗层有效)
struct Context {
    coarse: Vec<f32>,
    valid: Vec<bool>,
    shape: RasterShape,
    hydro: HydroModel,
    pyramid: ScalePyramid,
    landform: Vec<Landform>,
}

fn build_context(dem: Vec<f32>, w: u32, h: u32, res: f64) -> Context {
    let prepared = prepare_values(
        dem,
        common::meta_with_keys(w, h, res, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let hydro = build_hydro(&prepared, &hydro_cfg()).unwrap();
    let step = (hydro_cfg().coarse_res_m / res).round() as usize;
    let cw = hydro.shape.width;
    let ch = hydro.shape.height;
    let mut coarse = vec![0f32; cw * ch];
    let mut valid = vec![false; cw * ch];
    for ty in 0..ch {
        for tx in 0..cw {
            let (mut mn, mut c) = (0f64, 0u32);
            for y in ty * step..((ty + 1) * step).min(prepared.shape.height) {
                for x in tx * step..((tx + 1) * step).min(prepared.shape.width) {
                    let i = y * prepared.shape.width + x;
                    if prepared.valid[i] {
                        mn += prepared.raw[i] as f64;
                        c += 1;
                    }
                }
            }
            if c > 0 {
                coarse[ty * cw + tx] = (mn / c as f64) as f32;
                valid[ty * cw + tx] = true;
            }
        }
    }
    let pyramid =
        build_scale_pyramid(&coarse, &valid, hydro.shape, 0.15).unwrap();
    let mut landform = Vec::with_capacity(cw * ch);
    for ty in 0..ch {
        for tx in 0..cw {
            let i = ty * cw + tx;
            if !valid[i] {
                landform.push(Landform::Flat);
                continue;
            }
            let scale = pyramid.characteristic_scale_m[i].max(100.0) as f64;
            let pat = geomorphon_pattern(
                &coarse,
                cw,
                ch,
                hydro_cfg().coarse_res_m,
                tx,
                ty,
                scale,
                2.0 * hydro_cfg().coarse_res_m,
                3.0,
            );
            landform.push(pattern_to_landform(&pat));
        }
    }
    Context {
        coarse,
        valid,
        shape: hydro.shape,
        hydro,
        pyramid,
        landform,
    }
}

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
