//! 融合版坝子 v2(5m 原生判据, 固定路径 hhgq)
//! 运行: cargo run -p topo_core --example basin_fused --release

use topo_core::distance::edt_with_index;
use topo_core::filter::focal_relief;
use topo_core::geotiff::{self, GeoMeta};
use topo_core::hydro::fill_and_route;
use topo_core::pipeline::{dilate_round, erode_round, fill_interior_holes, fill_nodata, pct_of};
use topo_core::segment::label_connected;
use topo_core::terrain::{box_downsample, slope_horn_degrees};

const DEM_PATH: &str = r"G:\tif_features\county_feature\hhgq\dem.tif";
const OUT: &str = "target/basin_fused";

fn main() {
    std::fs::create_dir_all(OUT).unwrap();
    for acc_km2 in [0.3f64, 1.0, 2.0, 5.0, 10.0, 20.0] {
        run_with_acc(acc_km2);
    }
}

fn run_with_acc(acc_km2: f64) {
    let (dem5, meta) = geotiff::read_f32(DEM_PATH).unwrap();
    let (w5, h5) = (meta.width as usize, meta.height as usize);
    let res_f = meta.resolution();
    println!("DEM {w5}x{h5} @{res_f}m");
    let n5 = w5 * h5;

    // 河网(25m 提取, 上采样 5m)
    let (mut dem_c, cw, ch) = box_downsample(&dem5, w5, h5, res_f, 25.0);
    let res_c = 25.0;
    fill_nodata(&mut dem_c, cw, ch);
    let acc_th = (acc_km2 * 1e6 / (res_c * res_c)) as u32;
    let fr = fill_and_route(&dem_c, cw, ch, 99999.0);
    let river25: Vec<bool> = fr.acc.iter().map(|&a| a >= acc_th).collect();
    drop(fr);
    let scale = (res_c / res_f).round() as usize;
    let mut river5 = vec![false; n5];
    for y in 0..h5 {
        let sy = (y / scale).min(ch - 1);
        for x in 0..w5 {
            river5[y * w5 + x] = river25[sy * cw + (x / scale).min(cw - 1)];
        }
    }
    drop(river25);

    // 5m 层判据(对齐 legacy 语义: slope/relief 的 5x5 窗在 5m 上)
    let (src5, dist5) = edt_with_index(&river5, w5, h5);
    drop(river5);
    let slope5 = slope_horn_degrees(&dem5, w5, h5, res_f);
    let relief5 = focal_relief(&dem5, w5, h5, 5);
    let mut cand: Vec<bool> = vec![false; n5];
    for i in 0..n5 {
        cand[i] = dist5[i] as f64 * res_f <= 500.0
            && (dem5[i] - dem5[src5[i] as usize]).abs() < 5.0
            && slope5[i] < 5.0
            && relief5[i] < 5.0;
    }
    drop(slope5);
    drop(relief5);
    drop(dist5);
    drop(src5);
    println!("5m 候选 {:.2}%", 100.0 * cand.iter().filter(|b| **b).count() as f64 / n5 as f64);

    // 下采样 25m 做对象级处理
    let mut cand_c = vec![false; cw * ch];
    for y in 0..ch {
        for x in 0..cw {
            // 25m 像元内 5m 候选占比过半才算(抗噪)
            let y0 = y * scale;
            let x0 = x * scale;
            let mut cnt = 0;
            let mut tot = 0;
            for dy in 0..scale {
                for dx in 0..scale {
                    let yy = (y0 + dy).min(h5 - 1);
                    let xx = (x0 + dx).min(w5 - 1);
                    tot += 1;
                    if cand[yy * w5 + xx] {
                        cnt += 1;
                    }
                }
            }
            cand_c[y * cw + x] = cnt * 2 > tot;
        }
    }
    drop(cand);

    // 桥接闭运算 50m + 面积 5000m² + 内部起伏<=15m
    let r_b = 50.0 / res_c;
    dilate_round(&mut cand_c, cw, ch, r_b);
    erode_round(&mut cand_c, cw, ch, r_b);
    let (lab, n_obj) = label_connected(&cand_c, cw, ch);
    drop(cand_c);
    let mut vals: Vec<Vec<f32>> = vec![Vec::new(); n_obj + 1];
    for y in 0..ch {
        for x in 0..cw {
            let l = lab[y * cw + x] as usize;
            if l > 0 {
                vals[l].push(dem_c[y * cw + x]);
            }
        }
    }
    let min_cells = (5000.0 / (res_c * res_c)).ceil() as u64;
    let mut kept = vec![false; n_obj + 1];
    let (mut p_a, mut p_i) = (0u32, 0u32);
    for k in 1..=n_obj {
        if (vals[k].len() as u64) < min_cells {
            continue;
        }
        p_a += 1;
        let mut v = std::mem::take(&mut vals[k]);
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if pct_of(&v, 0.95) - pct_of(&v, 0.05) <= 15.0 {
            p_i += 1;
            kept[k] = true;
        }
    }
    println!("对象 {n_obj}, 面积达标 {p_a}, +内部起伏 {p_i}");
    let basin_c: Vec<bool> = lab.iter().map(|&l| l > 0 && kept[l as usize]).collect();
    drop(lab);
    drop(dem_c);

    // 上采样 5m + 填洞 + 平滑 50m
    let mut basin = vec![false; n5];
    for y in 0..h5 {
        let sy = (y / scale).min(ch - 1);
        for x in 0..w5 {
            basin[y * w5 + x] = basin_c[sy * cw + (x / scale).min(cw - 1)];
        }
    }
    drop(basin_c);
    fill_interior_holes(&mut basin, w5, h5);
    let r_sm = 50.0 / res_f;
    dilate_round(&mut basin, w5, h5, r_sm);
    erode_round(&mut basin, w5, h5, r_sm);

    let cnt = basin.iter().filter(|b| **b).count();
    let km2 = cnt as f64 * res_f * res_f / 1e6;
    println!("融合版坝子 {km2:.2} km² ({:.2}%)", 100.0 * km2 / (n5 as f64 * res_f * res_f / 1e6));
    let meta_out = GeoMeta::from_origin(w5 as u32, h5 as u32, 0.0, 0.0, res_f);
    let b8: Vec<u8> = basin.iter().map(|&b| b as u8).collect();
    let mut cmap = [[0u8; 3]; 256];
    cmap[1] = [51, 178, 229];
    geotiff::write_u8_cmap(format!("{OUT}/fused_{acc_km2}.tif"), &meta_out, &b8, &cmap).unwrap();
}
