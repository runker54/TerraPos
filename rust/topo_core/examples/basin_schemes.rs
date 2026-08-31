//! 坝子提取三方案对比(hhgq 真实数据, 固定路径)
//! 方案一 HAND 河谷低平带 | 方案二 填洼深度洼地 | 方案三 TPI 负地形
//! 运行: cargo run -p topo_core --example basin_schemes --release

use topo_core::distance::edt_with_index;
use topo_core::filter::{focal_mean, focal_relief};
use topo_core::geotiff::{self, GeoMeta};
use topo_core::hydro::fill_and_route;
use topo_core::pipeline::{
    closing_round, dilate_round, erode_round, fill_interior_holes, fill_nodata,
};
use topo_core::segment::label_connected;
use topo_core::terrain::{box_downsample, slope_horn_degrees};

const DEM_PATH: &str = r"G:\tif_features\county_feature\hhgq\dem.tif";
const OUT_BASE: &str = "target/basin_schemes";

/// 后处理+统计+落盘
fn finish(
    basin_c: &[bool],
    cw: usize,
    ch: usize,
    w5: usize,
    h5: usize,
    res_f: f64,
    scale: usize,
    meta_out: &GeoMeta,
    name: &str,
) {
    let mut b5 = vec![false; w5 * h5];
    for y in 0..h5 {
        let sy = (y / scale).min(ch - 1);
        for x in 0..w5 {
            b5[y * w5 + x] = basin_c[sy * cw + (x / scale).min(cw - 1)];
        }
    }
    let r_open = (30.0 / res_f).round();
    erode_round(&mut b5, w5, h5, r_open);
    dilate_round(&mut b5, w5, h5, r_open);
    let (lab, ln) = label_connected(&b5, w5, h5);
    let mut sizes = vec![0u64; ln + 1];
    for &l in lab.iter() {
        if l > 0 {
            sizes[l as usize] += 1;
        }
    }
    let min_cells = (150.0 * 666.67 / (res_f * res_f)).round() as u64; // 150 亩
    for i in 0..w5 * h5 {
        let l = lab[i];
        if l > 0 && sizes[l as usize] < min_cells {
            b5[i] = false;
        }
    }
    fill_interior_holes(&mut b5, w5, h5);
    let r_sm = (100.0 / res_f).round();
    dilate_round(&mut b5, w5, h5, r_sm);
    erode_round(&mut b5, w5, h5, r_sm);

    let n5 = b5.iter().filter(|b| **b).count();
    let (lab2, ln2) = label_connected(&b5, w5, h5);
    let mut sz2 = vec![0u64; ln2 + 1];
    for &l in lab2.iter() {
        if l > 0 {
            sz2[l as usize] += 1;
        }
    }
    let n_poly = (1..=ln2).filter(|&k| sz2[k] > 0).count();
    let km2 = n5 as f64 * res_f * res_f / 1e6;
    let pct = 100.0 * km2 / (w5 as f64 * h5 as f64 * res_f * res_f / 1e6);
    println!(
        "{:<16} {:>9} 像元 | {:>8.2} km² ({:.2}%) | 盆底 {} 个",
        name, n5, km2, pct, n_poly
    );
    let b8: Vec<u8> = b5.iter().map(|&b| b as u8).collect();
    let mut cmap = [[0u8; 3]; 256];
    cmap[1] = [51, 178, 229];
    let _ = geotiff::write_u8_cmap(format!("{OUT_BASE}/{name}.tif"), meta_out, &b8, &cmap);
}

fn main() {
    // ---------- 公共底图: 5m DEM + 25m 全填洼分析层 ----------
    std::fs::create_dir_all(OUT_BASE).unwrap();
    let (dem5, meta_raw) = geotiff::read_f32(DEM_PATH).unwrap();
    let (w5, h5) = (meta_raw.width as usize, meta_raw.height as usize);
    let res_f = meta_raw.resolution();
    let meta_out = GeoMeta::from_origin(w5 as u32, h5 as u32, 0.0, 0.0, res_f);
    println!("DEM {w5}x{h5} @{res_f}m");

    let (mut dem_c, cw, ch) = box_downsample(&dem5, w5, h5, res_f, 25.0);
    let res_c = 25.0;
    fill_nodata(&mut dem_c, cw, ch);
    // 全填洼(无 z-limit): 消除 DEM 异常值/洼地噪声, 三方案统一底图
    let fr_full = fill_and_route(&dem_c, cw, ch, 99999.0);
    let filled_c = fr_full.filled.clone();
    println!("25m 分析层 {cw}x{ch}, 全填洼完成");

    // 派生因子(全填洼面上计算, 消除异常值)
    let slope_c = slope_horn_degrees(&filled_c, cw, ch, res_c);
    let relief_c = focal_relief(&filled_c, cw, ch, 5);
    let mean2k = focal_mean(&filled_c, cw, ch, (2000.0 / res_c) as usize); // 2km TPI
    let acc_b = (0.3 * 1e6 / (res_c * res_c)) as u32; // 0.3km² 主干河网
    let river_b: Vec<bool> = fr_full.acc.iter().map(|&a| a >= acc_b).collect();
    let (rb_src, _) = edt_with_index(&river_b, cw, ch);
    let handb_c: Vec<f32> = filled_c
        .iter()
        .zip(rb_src.iter())
        .map(|(&z, &s)| (z - filled_c[s as usize]).abs())
        .collect();
    drop(rb_src);
    drop(river_b);
    drop(dem5);

    // ---------- 三方案候选(粗层判定) ----------
    let mut basin1 = vec![false; cw * ch]; // 方案一 HAND 低平带
    let mut basin2 = vec![false; cw * ch]; // 方案二 填洼深度
    let mut basin3 = vec![false; cw * ch]; // 方案三 TPI 负地形
    for i in 0..cw * ch {
        let depth = filled_c[i] - dem_c[i];
        let tpi = filled_c[i] - mean2k[i];
        basin1[i] = slope_c[i] < 4.0 && handb_c[i] < 10.0 && relief_c[i] < 4.0;
        basin2[i] = depth > 0.5 && slope_c[i] < 6.0 && relief_c[i] < 5.0;
        basin3[i] = tpi < -10.0 && slope_c[i] < 6.0 && relief_c[i] < 5.0;
    }
    drop(mean2k);
    drop(filled_c);
    drop(dem_c);

    // ---------- 通用后处理与输出 ----------
    let scale = (res_c / res_f).round() as usize;
    finish(&basin1, cw, ch, w5, h5, res_f, scale, &meta_out, "scheme1_hand(HAND低平带)");
    finish(&basin2, cw, ch, w5, h5, res_f, scale, &meta_out, "scheme2_depth(填洼深度)");
    finish(&basin3, cw, ch, w5, h5, res_f, scale, &meta_out, "scheme3_tpi(TPI负地形)");
    println!("完成, 栅格已写入 {OUT_BASE}/");
}
