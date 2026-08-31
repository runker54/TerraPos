//! 真实数据(hhgq 5m)全流程验证(方案三 TPI 负地形坝子)
//! 运行: cargo run -p topo_core --example run_hhgq --release

use std::sync::atomic::AtomicBool;
use topo_core::pipeline::{run, Params};

fn main() {
    let params = Params {
        dem_path: r"G:\tif_features\county_feature\hhgq\dem.tif".into(),
        out_dir: "target/hhgq_out".into(),
        ..Default::default()
    };
    let cancelled = AtomicBool::new(false);
    let t0 = std::time::Instant::now();
    let out = run(&params, &|p| {
        println!("[{:>5.1}%] {:<6} {}", p.pct, p.stage, p.msg);
        true
    }, &cancelled)
    .expect("pipeline failed");
    println!("耗时: {:.1}s", t0.elapsed().as_secs_f32());
    println!("{}", out.report);
}
