//! 参数效应契约: 每个暴露参数必须对其中间层或结果产生可观察的影响。

mod common;

use common::{full_chain, Context};
use std::sync::atomic::AtomicBool;
use topo_core::basin::{detect_basins, BasinConfig};
use topo_core::hydro::{build_hydro, HydroConfig};
use topo_core::input::{prepare_values, RasterShape};
use topo_core::pipeline::{
    run, BasinTendency, Params, PrecisionPreset,
};
use topo_core::ridge::build_ridges;
use topo_core::scale::build_scale_pyramid;
use topo_core::slope_unit::{build_slope_geometry, build_slope_units};

fn chain_with(
    dem: Vec<f32>,
    w: u32,
    h: u32,
    res: f64,
    hydro_z_limit: f32,
    growth: f32,
    low_relief: f32,
) -> (Context, topo_core::ridge::RidgeModel, topo_core::slope_unit::SlopeUnits) {
    let prepared = prepare_values(
        dem,
        common::meta_with_keys(w, h, res, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let hcfg = HydroConfig {
        coarse_res_m: 25.0,
        z_limit_m: hydro_z_limit,
        stream_areas_km2: [0.05, 0.20, 1.00, 5.00],
    };
    let hydro = build_hydro(&prepared, &hcfg).unwrap();
    let cw = hydro.shape.width;
    let ch = hydro.shape.height;
    let step = (hcfg.coarse_res_m / res).ceil() as usize;
    let mut coarse = vec![0f32; cw * ch];
    let mut valid = vec![false; cw * ch];
    for ty in 0..ch {
        let y0 = ty * step;
        let y1 = ((ty + 1) * step).min(prepared.shape.height);
        for tx in 0..cw {
            let x0 = tx * step;
            let x1 = ((tx + 1) * step).min(prepared.shape.width);
            let (mut sum, mut c) = (0f64, 0u32);
            for y in y0..y1 {
                for x in x0..x1 {
                    let i = y * prepared.shape.width + x;
                    if prepared.valid[i] {
                        sum += prepared.raw[i] as f64;
                        c += 1;
                    }
                }
            }
            if c > 0 {
                coarse[ty * cw + tx] = (sum / c as f64) as f32;
                valid[ty * cw + tx] = true;
            }
        }
    }
    let pyramid = build_scale_pyramid(&coarse, &valid, hydro.shape, growth).unwrap();
    let landform: Vec<topo_core::geomorphon::Landform> = (0..cw * ch)
        .map(|i| {
            if !valid[i] {
                return topo_core::geomorphon::Landform::Flat;
            }
            let sc = pyramid.characteristic_scale_m[i].max(100.0) as f64;
            let pat = topo_core::geomorphon::geomorphon_pattern(
                &coarse,
                cw,
                ch,
                25.0,
                i % cw,
                i / cw,
                sc,
                50.0,
                3.0,
            );
            topo_core::geomorphon::pattern_to_landform(&pat)
        })
        .collect();
    let ctx = Context {
        coarse,
        valid,
        shape: hydro.shape,
        hydro,
        pyramid,
        landform,
    };
    let ridges =
        build_ridges(&ctx.coarse, &ctx.valid, &ctx.hydro, &ctx.pyramid, &ctx.landform).unwrap();
    let units = build_slope_units(&ctx.coarse, &ctx.valid, &ctx.hydro, &ridges).unwrap();
    let _ = low_relief;
    (ctx, ridges, units)
}

/// 预设映射: Fast/Standard/Detailed 解析出不同的粗分辨率与尺度族
#[test]
fn precision_changes_coarse_resolution_and_scale_family() {
    let mut p = Params::default();
    p.precision = PrecisionPreset::Fast;
    let fast = p.resolved_with(5.0);
    assert_eq!(fast.coarse_res_m, 50.0, "Fast 应为 max(res,50)");
    assert_eq!(fast.scale_min_m, 250.0);
    let mut detailed = p.clone();
    detailed.precision = PrecisionPreset::Detailed;
    let d = detailed.resolved_with(5.0);
    assert_eq!(d.coarse_res_m, 10.0, "Detailed 应为 max(res,10)");
    assert_eq!(d.scale_max_m, 8000.0);
}

/// Strict->Loose 单调放宽盆地候选面积
#[test]
fn basin_tendency_monotone_relaxes_candidates() {
    let (dem, shape) = common::broad_basin(401, 401, 10.0);
    let (ctx, ridges, units, geom, morph) = full_chain(dem, shape.width as u32, shape.height as u32, shape.resolution_m);
    let cand = |t: BasinTendency| -> usize {
        detect_basins(
            &ctx.coarse,
            &ctx.valid,
            ctx.shape,
            &ctx.hydro,
            &units,
            &geom,
            &ctx.pyramid,
            &morph,
            &BasinConfig { tendency: t, min_area_m2: 66_666.67 },
        )
        .unwrap()
        .candidate
        .iter()
        .filter(|&&b| b)
        .count()
    };
    let strict = cand(BasinTendency::Strict);
    let loose = cand(BasinTendency::Loose);
    assert!(loose >= strict, "Loose 候选应 >= Strict: {loose} vs {strict}");
}

/// 增大 basin_min_area_m2 只减少接受盆地面积
#[test]
fn raising_min_area_monotone() {
    let (dem, shape) = common::broad_basin(401, 401, 10.0);
    let (ctx, ridges, units, geom, morph) = full_chain(dem, shape.width as u32, shape.height as u32, shape.resolution_m);
    let acc = |min_area: f64| -> usize {
        detect_basins(
            &ctx.coarse,
            &ctx.valid,
            ctx.shape,
            &ctx.hydro,
            &units,
            &geom,
            &ctx.pyramid,
            &morph,
            &BasinConfig { tendency: BasinTendency::Standard, min_area_m2: min_area },
        )
        .unwrap()
        .mask
        .iter()
        .filter(|&&b| b)
        .count()
    };
    let a1 = acc(66_666.67);
    let a2 = acc(5.0e7);
    assert!(a1 > 0, "应存在接受盆地");
    assert_eq!(a2, 0, "超阈值后应无接受盆地");
}

/// hydro_z_limit 改变深洼判别(更严的限深产生更多深洼)
#[test]
fn hydro_z_limit_changes_deep_sink_detection() {
    let deep = |z_limit: f32| -> usize {
        let (dem, _) = common::closed_pit(201, 201, 10.0);
        let prepared = prepare_values(
            dem,
            common::meta_with_keys(201, 201, 10.0, 1, 9001),
            &Default::default(),
        )
        .unwrap();
        let hcfg = HydroConfig {
            coarse_res_m: 25.0,
            z_limit_m: z_limit,
            stream_areas_km2: [0.05, 0.20, 1.00, 5.00],
        };
        build_hydro(&prepared, &hcfg)
            .unwrap()
            .deep_sink
            .iter()
            .filter(|&&b| b)
            .count()
    };
    let strict = deep(5.0);
    let lenient = deep(80.0);
    assert!(strict > lenient, "更严的限深应产生更多深洼: {strict} vs {lenient}");
}

/// 增长率阈值改变特征尺度(更高阈值 -> 更早收敛 -> 更小尺度)
#[test]
fn growth_threshold_changes_characteristic_scale() {
    let median_scale = |growth: f32| -> f32 {
        let (z, shape) = common::conical_hill(301, 301, 25.0);
        let valid = vec![true; z.len()];
        let p = build_scale_pyramid(&z, &valid, shape, growth).unwrap();
        let mut v: Vec<f32> = p
            .characteristic_scale_m
            .iter()
            .copied()
            .filter(|s| *s > 0.0)
            .collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    let tight = median_scale(0.05);
    let loose = median_scale(0.40);
    assert!(
        loose <= tight,
        "更高增长率阈值不应给出更大的特征尺度: {loose} vs {tight}"
    );
}

/// low_relief_m 改变相对位置的融合权重(低起伏判据变化 -> q 分布变化)
#[test]
fn low_relief_threshold_changes_relative_position() {
    let mean_q = |low_relief: f32| -> f32 {
        let (dem, shape) = common::v_valley(201, 201, 10.0);
        let (ctx, _r, _u) = chain_with(dem, shape.width as u32, shape.height as u32, shape.resolution_m, 15.0, 0.15, low_relief);
        let geom = build_slope_geometry(
            &ctx.coarse,
            &ctx.valid,
            &units_helper(&ctx),
            &ctx.hydro,
            &ctx.pyramid,
            low_relief,
        )
        .unwrap();
        let s: f64 = geom
            .relative_position
            .iter()
            .copied()
            .filter(|v| v.is_finite())
            .map(|v| v as f64)
            .sum();
        let cnt = geom.relative_position.iter().filter(|v| v.is_finite()).count();
        if cnt == 0 { 0.0 } else { (s / cnt as f64) as f32 }
    };
    let a = mean_q(5.0);
    let b = mean_q(200.0);
    assert!(
        (a - b).abs() > 1e-4,
        "low_relief_m 未改变相对位置分布: {a} vs {b}"
    );
}

fn units_helper(ctx: &Context) -> topo_core::slope_unit::SlopeUnits {
    let ridges =
        build_ridges(&ctx.coarse, &ctx.valid, &ctx.hydro, &ctx.pyramid, &ctx.landform).unwrap();
    build_slope_units(&ctx.coarse, &ctx.valid, &ctx.hydro, &ridges).unwrap()
}

/// hill_elevation_max_m 决定丘陵/山地归类边界
#[test]
fn hill_elevation_threshold_switches_classes() {
    let n = 4;
    let elev = vec![400.0f32, 400.0, 600.0, 600.0];
    let valid = vec![true; n];
    let basin = topo_core::basin::BasinResult::default();
    let positions = topo_core::slope_position::Memberships {
        upper: vec![0.0; n],
        middle: vec![0.0; n],
        lower: vec![1.0; n],
        raw: vec![topo_core::slope_position::SlopePosition::Lower; n],
        confidence: vec![1.0; n],
    };
    let fc = topo_core::postprocess::compose_codes(&elev, &valid, &basin, &positions, 500.0).unwrap();
    assert_eq!(fc.terrain[0], 5, "400m 应为丘陵下部");
    assert_eq!(fc.terrain[2], 8, "600m 应为山地坡下");
    let fc2 = topo_core::postprocess::compose_codes(&elev, &valid, &basin, &positions, 300.0).unwrap();
    assert_eq!(fc2.terrain[0], 8, "300m 阈值下 400m 应为山地");
}

/// 嵌套河网面积等级改变 stream_level 分布
#[test]
fn stream_area_thresholds_change_stream_levels() {
    let level_count = |areas: [f64; 4]| -> [usize; 5] {
        let (dem, _) = common::v_valley(401, 401, 10.0);
        let prepared = prepare_values(
            dem,
            common::meta_with_keys(401, 401, 10.0, 1, 9001),
            &Default::default(),
        )
        .unwrap();
        let hcfg = HydroConfig {
            coarse_res_m: 25.0,
            z_limit_m: 15.0,
            stream_areas_km2: areas,
        };
        let hydro = build_hydro(&prepared, &hcfg).unwrap();
        let mut cnt = [0usize; 5];
        for &lv in &hydro.stream_level {
            cnt[lv as usize] += 1;
        }
        cnt
    };
    let fine = level_count([0.05, 0.20, 1.00, 5.00]);
    let coarse = level_count([0.50, 2.00, 10.0, 50.0]);
    let fine_total: usize = fine[1..].iter().sum();
    let coarse_total: usize = coarse[1..].iter().sum();
    assert!(
        fine_total > coarse_total,
        "更细的最低阈值应产生更多河网像元: {:?} vs {:?}",
        fine,
        coarse
    );
}

/// write_diagnostics 只影响文件产出, 不改变编码结果
#[test]
fn write_diagnostics_changes_files_only() {
    let dir = std::env::temp_dir().join("topo_param_effect");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    // 写一个小的投影米制 DEM
    let (dem, shape) = common::nested_ridges(257, 257, 10.0);
    let p = dir.join("dem.tif");
    let meta = common::meta_with_keys(shape.width as u32, shape.height as u32, shape.resolution_m, 1, 9001);
    topo_core::geotiff::write_f32(&p, &meta, &dem).unwrap();

    let make_params = |diag: bool, out: String| -> Params {
        Params {
            dem_path: p.to_string_lossy().into_owned(),
            out_dir: out,
            precision: PrecisionPreset::Standard,
            basin_tendency: BasinTendency::Standard,
            basin_min_area_m2: 66_666.67,
            postprocess_strength: 1.0,
            write_diagnostics: diag,
            advanced: Default::default(),
        }
    };
    let cancelled = AtomicBool::new(false);
    eprintln!("STEP before o1");
    let o1 = run(
        &make_params(true, dir.join("with_diag").to_string_lossy().into_owned()),
        &|_| true,
        &cancelled,
    );
    eprintln!("STEP o1 done: {:?}", o1.as_ref().map(|o| o.terrain.len()));
    let o1 = o1.unwrap();
    let o2 = run(
        &make_params(false, dir.join("no_diag").to_string_lossy().into_owned()),
        &|_| true,
        &cancelled,
    )
    .unwrap();
    assert_eq!(o1.terrain, o2.terrain, "诊断开关不得改变编码");
    assert!(
        dir.join("with_diag").join("diagnostics").exists(),
        "开启诊断应产出 diagnostics 目录"
    );
    assert!(
        !dir.join("no_diag").join("diagnostics").exists(),
        "关闭诊断不应产出 diagnostics 目录"
    );
}
