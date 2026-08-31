//! 端到端运行: DEM 全流程 + 耗时统计
//! 用法: run_full [dem_path] [out_dir]   (缺省跑内置大 DEM)
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;
use topo_core::pipeline::{run, Params, Progress};

fn main() {
    let mut args = std::env::args().skip(1);
    let dem_path = args
        .next()
        .unwrap_or_else(|| r"E:\zcode_worker\Topographic\data\dem.tif".into());
    let out_dir = args
        .next()
        .unwrap_or_else(|| r"E:\zcode_worker\Topographic\app\app_out".into());
    let params = Params { dem_path, out_dir, ..Params::default() };

    let cancel = AtomicBool::new(false);
    let t0 = Instant::now();
    let result = run(&params, &|p: Progress| {
        println!("[{:5.1}%] {}: {}", p.pct, p.stage, p.msg);
        true
    }, &cancel);
    match result {
        Ok(out) => {
            println!("\n=== 完成, 总耗时 {:?} ===", t0.elapsed());
            println!("stats = {:?}", out.stats);
            println!("{}", out.report);
        }
        Err(e) => {
            eprintln!("失败: {:?} (耗时 {:?})", e, t0.elapsed());
            std::process::exit(1);
        }
    }
    let _ = cancel.load(Ordering::Relaxed);
}
