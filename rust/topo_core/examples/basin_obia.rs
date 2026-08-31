//! 山间盆地识别 v4: 面向对象(OBIA)复合判据 —— hhgq 实测
//! 平坦面对象分割 -> 对象级判据(面积/内部起伏/包围度/内切宽度) -> 后处理
//! 运行: cargo run -p topo_core --example basin_obia --release

use topo_core::distance::edt_with_index;
use topo_core::filter::{focal_relief};
use topo_core::geotiff::{self, GeoMeta};
use topo_core::hydro::fill_and_route;
use topo_core::pipeline::{dilate_round, erode_round, fill_interior_holes, fill_nodata};
use topo_core::segment::label_connected;
use topo_core::terrain::{box_downsample, slope_horn_degrees};

const DEM_PATH: &str = r"G:\tif_features\county_feature\hhgq\dem.tif";
const OUT: &str = "target/basin_obia";

/// 分位数(升序 vec, p∈[0,1])
fn pct(sorted: &[f32], p: f64) -> f32 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn main() {
    std::fs::create_dir_all(OUT).unwrap();
    let (dem5, meta) = geotiff::read_f32(DEM_PATH).unwrap();
    let (w5, h5) = (meta.width as usize, meta.height as usize);
    let res_f = meta.resolution();
    println!("DEM {w5}x{h5} @{res_f}m");

    // 25m 粗层 + 全填洼
    let (mut dem_c, cw, ch) = box_downsample(&dem5, w5, h5, res_f, 25.0);
    let res_c = 25.0;
    fill_nodata(&mut dem_c, cw, ch);
    let filled = fill_and_route(&dem_c, cw, ch, 99999.0).filled;
    println!("粗层 {cw}x{ch} 全填洼完成");

    // 因子
    let slope = slope_horn_degrees(&filled, cw, ch, res_c);
    let relief = focal_relief(&filled, cw, ch, 5);

    // ---- 1. 平坦候选 + 桥接闭运算(100m, legacy buffer 语义: 沟渠/道路碎片重组) ----
    let mut flat: Vec<bool> = (0..cw * ch)
        .map(|i| slope[i] < 6.0 && relief[i] < 5.0)
        .collect();
    dilate_round(&mut flat, cw, ch, 4.0); // 100m
    erode_round(&mut flat, cw, ch, 4.0);
    let n_flat = flat.iter().filter(|b| **b).count();
    println!("平坦候选 {:.1}%", 100.0 * n_flat as f64 / (cw * ch) as f64);

    // ---- 2. 对象分割 ----
    let (lab, n_obj) = label_connected(&flat, cw, ch);
    drop(flat);
    println!("平坦面对象 {n_obj} 个");

    // ---- 3. 对象级判据 ----
    // 参数(hhgq 校准)
    let min_cells = (150.0 * 666.67 / (res_c * res_c)) as u64; // 150 亩 = 160 px
    let th_inner: f32 = 15.0; // 内部起伏 P95-P5 <= 15m
    let th_surround: f32 = 10.0; // 环带(200m) P75 - 对象 P50 >= 10m
    let th_md_px: f32 = (250.0 / 2.0 / res_c) as f32; // 内切全宽>=250m -> 半径>=5px

    // 内切半径场: 到背景(补集)的距离才是内切半径(EDT 源必须取补集)
    let flat_mask: Vec<bool> = lab.iter().map(|&l| l > 0).collect();
    let inv: Vec<bool> = flat_mask.iter().map(|b| !b).collect();
    let (_, md) = edt_with_index(&inv, cw, ch);
    drop(flat_mask);
    drop(inv);

    // 对象像元/环带高程收集(一次遍历)
    let mut obj_vals: Vec<Vec<f32>> = vec![Vec::new(); n_obj + 1];
    let mut obj_cells: Vec<u64> = vec![0; n_obj + 1];
    let mut obj_md: Vec<f32> = vec![0.0; n_obj + 1];
    for i in 0..cw * ch {
        let l = lab[i] as usize;
        if l > 0 {
            obj_vals[l].push(filled[i]);
            obj_cells[l] += 1;
            if md[i] > obj_md[l] {
                obj_md[l] = md[i];
            }
        }
    }
    drop(md);

    // 环带: 对象并集膨胀 8px(200m) 减对象
    let obj_mask: Vec<bool> = lab.iter().map(|&l| l > 0).collect();
    let cands: Vec<usize> = (1..=n_obj)
        .filter(|&k| obj_cells[k] >= min_cells)
        .collect();
    println!("面积达标(>=150亩)对象 {} 个, 计算包围度...", cands.len());
    let mut kept = vec![false; n_obj + 1];
    let mut pass_cnt = [0u32; 4]; // [面积, 起伏, 包围, 宽度] 逐级通过数
    pass_cnt[0] = cands.len() as u32;
    for &k in &cands {
        let mut v = std::mem::take(&mut obj_vals[k]);
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let inner = pct(&v, 0.95) - pct(&v, 0.05);
        let med = pct(&v, 0.5);
        obj_vals[k] = v;
        if inner > th_inner {
            continue;
        }
        pass_cnt[1] += 1;
        // 该对象自己的环带
        let m: Vec<bool> = lab.iter().map(|&l| l == k as i32).collect();
        let mut r = m.clone();
        dilate_round(&mut r, cw, ch, 8.0);
        let mut rv = Vec::new();
        for i in 0..cw * ch {
            if r[i] && !m[i] {
                rv.push(filled[i]);
            }
        }
        if rv.is_empty() {
            continue;
        }
        rv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let surround = pct(&rv, 0.75) - med;
        if surround < th_surround {
            continue;
        }
        pass_cnt[2] += 1;
        if obj_md[k] >= th_md_px {
            pass_cnt[3] += 1;
            kept[k] = true;
        }
    }

    // 诊断: 面积达标对象的内切宽分布
    {
        let mut ws: Vec<f64> = cands.iter().map(|&k| obj_md[k] as f64 * 2.0 * res_c).collect();
        ws.sort_by(|a: &f64, b: &f64| a.partial_cmp(b).unwrap());
        let n = ws.len();
        if n > 0 {
            println!(
                "对象内切全宽: 中位 {:.0}m, p25 {:.0}m, p75 {:.0}m, max {:.0}m",
                ws[n / 2],
                ws[n / 4],
                ws[3 * n / 4],
                ws[n - 1]
            );
        }
    }

    // ---- 4. 合格对象并集 ----
    let mut basin_c: Vec<bool> = lab.iter().map(|&l| l > 0 && kept[l as usize]).collect();
    drop(lab);

    // ---- 5. 上采样 5m + 填洞 + 平滑 ----
    let scale = (res_c / res_f).round() as usize;
    let mut basin = vec![false; w5 * h5];
    for y in 0..h5 {
        let sy = (y / scale).min(ch - 1);
        for x in 0..w5 {
            basin[y * w5 + x] = basin_c[sy * cw + (x / scale).min(cw - 1)];
        }
    }
    drop(basin_c);
    fill_interior_holes(&mut basin, w5, h5);
    let r_sm = (100.0 / res_f).round();
    dilate_round(&mut basin, w5, h5, r_sm);
    erode_round(&mut basin, w5, h5, r_sm);

    // ---- 6. 统计与落盘 ----
    let n5 = basin.iter().filter(|b| **b).count();
    let km2 = n5 as f64 * res_f * res_f / 1e6;
    let (_, ln2) = label_connected(&basin, w5, h5);
    println!(
        "判据瀑布: 面积达标 {} -> +内部起伏<={}m {} -> +包围>={}m {} -> +内切宽>={}m {}",
        pass_cnt[0], th_inner, pass_cnt[1], th_surround, pass_cnt[2], 250, pass_cnt[3]
    );
    println!("坝子 {km2:.2} km² ({:.2}%), 成品层连通域 {ln2} 个",
        100.0 * km2 / (w5 as f64 * h5 as f64 * res_f * res_f / 1e6));
    let meta_out = GeoMeta::from_origin(w5 as u32, h5 as u32, 0.0, 0.0, res_f);
    let b8: Vec<u8> = basin.iter().map(|&b| b as u8).collect();
    let mut cmap = [[0u8; 3]; 256];
    cmap[1] = [51, 178, 229];
    geotiff::write_u8_cmap(format!("{OUT}/basin_obia.tif"), &meta_out, &b8, &cmap).unwrap();
    println!("已写出 {OUT}/basin_obia.tif");
}
