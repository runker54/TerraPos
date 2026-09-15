//! 批量验证 runner(无命令行参数, 路径常量在文件顶部):
//! 1. 扫描 data/synthetic/generated/*.tif 逐场景跑完整管线;
//! 2. 对 wide_valley 跑 base/±20% 参数扰动组;
//! 3. 跑 data/dem.tif 真 DEM。
//! 输出到 data/validation_out/, 供 python/validate_terrain_result.py 判定。

use std::sync::atomic::AtomicBool;
use topo_core::pipeline::{run, BasinTendency, Params, PrecisionPreset};

const SCENE_DIR: &str = "../data/synthetic/generated";
const VALIDATION_OUT: &str = "../data/validation_out";
const REAL_DEM: &str = "../data/dem.tif";

fn scene_params(dem: &str, out: &str, area_scale: f64, strength_scale: f64, z_scale: f64) -> Params {
    let mut adv = topo_core::pipeline::AdvancedParams::default();
    adv.hydro_z_limit_m = (15.0 * z_scale) as f32;
    Params {
        dem_path: dem.to_string(),
        out_dir: out.to_string(),
        precision: PrecisionPreset::Standard,
        basin_tendency: BasinTendency::Standard,
        basin_min_area_m2: 66_666.67 * area_scale,
        postprocess_strength: strength_scale,
        write_diagnostics: true,
        advanced: adv,
    }
}

fn run_one(params: &Params) -> Result<(), String> {
    let cancelled = AtomicBool::new(false);
    run(params, &|p| {
        print!("\r[{:>5.1}%] {}", p.pct, p.msg);
        true
    }, &cancelled)
    .map(|_| ())
    .map_err(|e| e.to_string())
}

fn main() {
    let mut failed = 0usize;
    // 1) 合成场景
    let scenes: Vec<_> = std::fs::read_dir(SCENE_DIR)
        .expect("场景目录不存在, 先运行 python/build_synthetic_dem.py")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension().map(|e| e == "tif").unwrap_or(false)
                && !p.file_name().unwrap().to_string_lossy().ends_with("_expected.tif")
                && !p.file_name().unwrap().to_string_lossy().starts_with("ref_")
        })
        .collect();
    for scene in &scenes {
        let name = scene.file_stem().unwrap().to_string_lossy().into_owned();
        let out_dir = format!("{VALIDATION_OUT}/synthetic/{name}");
        println!("\n== 场景 {name} ==");
        let params = scene_params(
            scene.to_string_lossy().as_ref(),
            &out_dir,
            1.0,
            1.0,
            1.0,
        );
        if let Err(e) = run_one(&params) {
            eprintln!("场景 {name} 失败: {e}");
            failed += 1;
        }
    }

    // 2) ±20% 参数扰动组(wide_valley)
    let wide = format!("{SCENE_DIR}/v_valley.tif");
    if std::path::Path::new(&wide).exists() {
        for (name, a, s, z) in [
            ("base", 1.0, 1.0, 1.0),
            ("minus", 0.8, 0.8, 0.8),
            ("plus", 1.2, 1.2, 1.2),
        ] {
            let out_dir = format!("{VALIDATION_OUT}/param/wide_valley/{name}");
            println!("\n== 参数组 {name} ==");
            let params = scene_params(&wide, &out_dir, a, s, z);
            if let Err(e) = run_one(&params) {
                eprintln!("参数组 {name} 失败: {e}");
                failed += 1;
            }
        }
    }

    // 3) 真 DEM
    if std::path::Path::new(REAL_DEM).exists() {
        println!("\n== 真 DEM ==");
        let params = Params {
            dem_path: REAL_DEM.to_string(),
            out_dir: format!("{VALIDATION_OUT}/dem"),
            ..Default::default()
        };
        if let Err(e) = run_one(&params) {
            eprintln!("真 DEM 失败: {e}");
            failed += 1;
        }
    }

    if failed > 0 {
        eprintln!("\n{failed} 个运行失败");
        std::process::exit(1);
    }
    println!("\n全部运行完成");
}
