//! 诊断: 打印合成挑战 DEM 的 prominence 种子及突出度
//! 运行: cargo run -p topo_core --example debug_prom --release

use topo_core::filter::focal_mean;
use topo_core::segment::prominence_map;

fn main() {
    // 与 compare_seeds 相同的合成 DEM
    let (w, h) = (600usize, 600usize);
    let hills = [
        (440.0, 410.0, 110.0, 5.0f64),
        (454.0, 422.0, 100.0, 5.0),
        (468.0, 410.0, 115.0, 6.0),
        (482.0, 422.0, 105.0, 5.0),
        (496.0, 410.0, 112.0, 6.0),
        (510.0, 422.0, 120.0, 6.0),
    ];
    let mut dem = vec![0f32; w * h];
    let gauss = |x: f64, y: f64, cx: f64, cy: f64, a: f64, s: f64| {
        let d2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
        a * (-d2 / (2.0 * s * s)).exp()
    };
    for y in 0..h {
        for x in 0..w {
            let (fx, fy) = (x as f64, y as f64);
            let mut z = 300.0 + fx * 0.01;
            z += gauss(fx, fy, 250.0, 250.0, 620.0, 50.0);
            z += gauss(fx, fy, 314.0, 314.0, 180.0, 12.0);
            z += gauss(fx, fy, 110.0, 120.0, 1200.0, 28.0);
            z += gauss(fx, fy, 124.0, 136.0, 1000.0, 25.0);
            for &(cx, cy, a, s) in &hills {
                z += gauss(fx, fy, cx, cy, a, s);
            }
            z += gauss(fx, fy, 170.0, 420.0, 30.0, 35.0);
            dem[y * w + x] = z as f32;
        }
    }
    let smooth = focal_mean(&dem, w, h, 5);
    let (prom, is_peak) = prominence_map(&smooth, w, h);
    let mut peaks: Vec<(f32, usize)> = (0..w * h)
        .filter(|&i| is_peak[i])
        .map(|i| (prom[i], i))
        .collect();
    peaks.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    println!("峰总数 {} (突出度降序, 前 25):", peaks.len());
    for &(p, i) in peaks.iter().take(25) {
        println!(
            "  ({:>3},{:>3}) 高程 {:>7.1} 突出度 {:>7.2}",
            i % w,
            i / w,
            smooth[i],
            p
        );
    }
}
