//! 端到端冒烟验证: 对样例 DEM 跑完整管线并打印分类统计
//! 运行: cargo run -p topo_core --example run_sample --release

use std::sync::atomic::AtomicBool;
use topo_core::pipeline::{run, Params};

fn main() {
    let params = Params {
        dem_path: "sample/sample_dem.tif".into(),
        out_dir: "target/verify_out".into(),
        ..Default::default()
    };
    let cancelled = AtomicBool::new(false);
    let t0 = std::time::Instant::now();
    let out = run(
        &params,
        &|p| {
            println!("[{:>5.1}%] {:<6} {}", p.pct, p.stage, p.msg);
            true
        },
        &cancelled,
    )
    .expect("pipeline failed");
    println!("\n耗时: {:.1}s", t0.elapsed().as_secs_f32());
    println!("{}", out.report);
    let mut cnt = [0u64; 256];
    for &c in &out.terrain {
        cnt[c as usize] += 1;
    }
    let total = out.terrain.len();
    println!("\n8类分布 (共 {} 像元):", total);
    for (code, name) in [
        (1u8, "山间盆地"),
        (2, "宽谷盆地"),
        (3, "丘陵下"),
        (4, "丘陵中"),
        (5, "丘陵上"),
        (6, "山地坡下"),
        (7, "山地坡中"),
        (8, "山地坡上"),
    ] {
        println!("  {} {:<6} {:>10}  {:>5.1}%", code, name, cnt[code as usize], 100.0 * cnt[code as usize] as f64 / total as f64);
    }
}
