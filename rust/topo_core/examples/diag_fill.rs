//! 填洼行为诊断(hhgq): 填充深度分布 + 假湖区与坝子的重叠
//! 运行: cargo run -p topo_core --example diag_fill --release

use topo_core::geotiff;
use topo_core::hydro::fill_and_route;
use topo_core::pipeline::fill_nodata;
use topo_core::terrain::box_downsample;

fn main() {
    let (dem5, meta) = geotiff::read_f32(r"G:\tif_features\county_feature\hhgq\dem.tif").unwrap();
    let (w5, h5) = (meta.width as usize, meta.height as usize);
    let res_f = meta.resolution();
    let (mut dem_c, cw, ch) = box_downsample(&dem5, w5, h5, res_f, 25.0);
    fill_nodata(&mut dem_c, cw, ch);

    // 无限制填洼(当前主管线行为)
    let fr = fill_and_route(&dem_c, cw, ch, 99999.0);
    let filled = fr.filled;

    // 填充深度 = filled - dem
    let n = cw * ch;
    let mut depths = Vec::with_capacity(n);
    for i in 0..n {
        let d = filled[i] - dem_c[i];
        if d > 0.01 {
            depths.push(d);
        }
    }
    depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let total = n as f64;
    let filled_pct = 100.0 * depths.len() as f64 / total;
    println!("粗层 {cw}x{ch} ({n} 像元)");
    println!("被填充像元: {:.1}% ({})", filled_pct, depths.len());
    if !depths.is_empty() {
        let q = |p: f64| depths[((depths.len() - 1) as f64 * p) as usize];
        println!("填充深度: 中位 {:.1}m, P90 {:.1}m, P99 {:.1}m, max {:.1}m",
            q(0.5), q(0.9), q(0.99), depths[depths.len()-1]);
        for th in [5.0f32, 10.0, 20.0, 50.0] {
            let c = depths.iter().filter(|&&d| d > th).count();
            println!("  深度>{}m: {} 像元 ({:.1}%) = {:.1} km²",
                th, c, 100.0 * c as f64 / total, c as f64 * 625.0 / 1e6);
        }
    }

    // 假湖判定: 填充深度 > 10m 的连通区, 取最大的几个看面积与位置
    let deep_fill: Vec<bool> = (0..n).map(|i| filled[i] - dem_c[i] > 10.0).collect();
    let (lab, ln) = topo_core::segment::label_connected(&deep_fill, cw, ch);
    let mut sizes = vec![0u64; ln + 1];
    for &l in lab.iter() {
        if l > 0 {
            sizes[l as usize] += 1;
        }
    }
    let mut order: Vec<(u64, usize)> = sizes.iter().enumerate().filter(|(k, _)| *k > 0)
        .map(|(k, &s)| (s, k)).collect();
    order.sort_by(|a, b| b.0.cmp(&a.0));
    println!("\n深填充(>10m)连通区 top5:");
    for (sz, k) in order.iter().take(5) {
        // 该连通区内的填充深度范围与原始地形起伏
        let mut dmax = 0f32;
        let mut rng = 0f32;
        let mut zmin = f32::MAX;
        let mut zmax = f32::MIN;
        for i in 0..n {
            if lab[i] as usize == *k {
                let d = filled[i] - dem_c[i];
                dmax = dmax.max(d);
                zmin = zmin.min(dem_c[i]);
                zmax = zmax.max(dem_c[i]);
                rng = rng.max(dmax);
            }
        }
        println!("  面积 {:.2} km²: 填深 max {:.0}m, 原始地面高差 {:.0}m ({}-{})",
            sz * 625 / 1_000_000, dmax, zmax - zmin, zmin as i32, zmax as i32);
    }

    // 与已识别坝子的重叠
    let (b5, _) = geotiff::read_f32("target/hhgq_out/basin_mask.tif").unwrap();
    let scale = (res_f / 25.0).round() as usize;
    let mut hit = 0u64;
    let mut tot = 0u64;
    for y in 0..h5 {
        let sy = (y / scale).min(ch - 1);
        for x in 0..w5 {
            let sx = (x / scale).min(cw - 1);
            if b5[y * w5 + x] > 0.0 {
                tot += 1;
                if filled[sy * cw + sx] - dem_c[sy * cw + sx] > 10.0 {
                    hit += 1;
                }
            }
        }
    }
    println!("\n已识别坝子中落在深填充区(>10m)的比例: {:.1}% ({} / {})",
        100.0 * hit as f64 / tot.max(1) as f64, hit, tot);
}
