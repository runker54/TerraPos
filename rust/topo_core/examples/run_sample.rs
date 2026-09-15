//! 端到端冒烟门: 对样例 DEM 跑完整管线, 校验分类契约(不使用目标面积比例)
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
    let valid = total - cnt[0] as usize; // NoData(0) 不进入分母

    // 分类契约断言(冒烟门): 保留码缺席 / 上中下齐备 / 无单类主导
    assert_eq!(cnt[2], 0, "代码 2 为保留码, 任何成果不得产生");
    assert!(cnt[3] + cnt[6] > 0, "缺少上部坡位(3/6)");
    assert!(cnt[4] + cnt[7] > 0, "缺少中部坡位(4/7)");
    assert!(cnt[5] + cnt[8] > 0, "缺少下部坡位(5/8)");
    for code in [3u8, 4, 5, 6, 7, 8] {
        let pct = 100.0 * cnt[code as usize] as f64 / valid as f64;
        assert!(
            pct <= 90.0,
            "非盆地类别 {code} 占有效像元 {pct:.1}% > 90%"
        );
    }

    println!("\n类别分布 (共 {total} 像元, 有效 {valid}):");
    for (code, name) in [
        (1u8, "山间/宽谷盆地"),
        (3, "丘陵上部"),
        (4, "丘陵中部"),
        (5, "丘陵下部"),
        (6, "山地坡上"),
        (7, "山地坡中"),
        (8, "山地坡下"),
    ] {
        println!(
            "  {} {:<10} {:>10}  {:>5.1}%",
            code,
            name,
            cnt[code as usize],
            100.0 * cnt[code as usize] as f64 / valid as f64
        );
    }
}
