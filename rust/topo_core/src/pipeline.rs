//! 参数模型 + 全流程编排（含进度回调与取消）

use crate::distance::edt_with_index;
use crate::error::{CoreError, Result};
use crate::filter::{focal_mean, focal_relief};
use crate::geotiff::{self, GeoMeta};
use crate::hydro::fill_and_route;
use crate::segment::label_connected;
use crate::terrain::{box_downsample, slope_horn_degrees};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// 山体单元种子提取模式
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default, Serialize, Deserialize)]
pub enum SeedMode {
    /// 单一固定窗口(峰顶间距 = window_peak_m)
    #[default]
    Fixed,
    /// 分亚类动态窗口(峰顶间距随地貌亚类自适应: 低丘窄、中山宽)
    Zoned,
    /// 地形突起度(峰顶 prominence >= seed_prominence_m, 无窗口参数)
    Prominence,
    /// 多尺度 TPI 特征尺度投票(尺度空间语义, 无种子; 正 TPI 连通域为个体)
    ScaleVote,
    /// 混合: 突出度语义 ∪ 距离语义(推荐)
    Hybrid,
}

/// 全部数值型指标参数（UI 表单一一对应）
#[derive(Debug, Clone)]
pub struct Params {
    pub dem_path: String,
    pub out_dir: String,
    /// 中间(大尺度)分析层分辨率(米), 默认 25; 需 >= 成品分辨率
    pub coarse_res: f64,
    /// 区域高程百分位上限(0-100): 低于邻域均值的幅度阈值
    pub basin_pctl_max: f64,
    /// 区域高程百分位邻域窗口(米)
    pub basin_pctl_win_m: f64,
    /// 坝子坡度上限(度)
    pub basin_slope_th: f64,
    /// 局部起伏(5x5@成品分辨率)上限(米)
    pub basin_relief_m: f64,
    /// 坝子最小保留面积(m²)
    pub basin_min_area_m2: f64,
    /// 对象内部起伏上限(米, P95-P5)
    pub basin_inner_relief_m: f64,
    /// 碎片桥接闭运算半径(米)
    pub basin_bridge_m: f64,
    /// 坝子内碎斑归并: 桥接域半径(米, 0=不归并)
    pub basin_merge_m: f64,
    /// 坝子内碎斑归并: 碎斑面积上限(m²)
    pub basin_merge_max_m2: f64,
    /// 坝子平滑距离(米)
    pub basin_smooth_m: f64,
    /// 坡位 geomorphon 线视搜索半径(米, 对齐 r.geomorphon search)
    pub slope_search_m: f64,
    /// 坡位线视起始跳过距离(米, 滤微起伏)
    pub slope_skip_m: f64,
    /// 平坡/中坡坡度分界(度, 对齐脚本阈值 5)
    pub slope_flat_deg: f64,
    /// 坡位小斑蚕食阈值(m², 对齐脚本 200 像元@5m)
    pub slope_min_patch_m2: f64,
    /// 丘陵海拔上限(米)
    pub hill_z_max: f64,
    /// 丘陵亚类起伏度窗口(米)
    pub relief_subclass_win: f64,
    /// 低丘起伏度上限(米)
    pub relief_low_hill: f64,
    /// 众数滤波轮数
    pub mode_filter_iter: usize,
    /// 最小图斑(m², 成品层)
    pub min_patch_m2: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            dem_path: String::new(),
            out_dir: String::new(),
            coarse_res: 25.0,
            basin_pctl_max: 30.0,
            basin_pctl_win_m: 2000.0,
            basin_slope_th: 5.0,
            basin_relief_m: 5.0,
            basin_min_area_m2: 66_666.67, // 100 亩
            basin_inner_relief_m: 30.0,
            basin_bridge_m: 50.0,
            basin_merge_m: 100.0,
            basin_merge_max_m2: 20_000.0,
            basin_smooth_m: 50.0,
            slope_search_m: 1000.0,
            slope_skip_m: 50.0,
            slope_flat_deg: 5.0,
            slope_min_patch_m2: 5000.0,
            hill_z_max: 500.0,
            relief_subclass_win: 2000.0,
            relief_low_hill: 200.0,
            mode_filter_iter: 1,
            min_patch_m2: 10_000.0,
        }
    }
}

pub struct Progress {
    pub stage: String,
    pub pct: f32,
    pub msg: String,
}

pub struct Outputs {
    pub terrain: Vec<u8>,
    pub subclass: Vec<u8>,
    pub meta5: GeoMeta,
    pub report: String,
    /// (类编码, 面积km²), 按业务顺序
    pub stats: Vec<(u8, f64)>,
}

/// 最近有效值填充 nodata
pub fn fill_nodata(dem: &mut [f32], w: usize, h: usize) -> bool {
    let invalid: Vec<bool> = dem.iter().map(|v| !v.is_finite()).collect();
    if !invalid.iter().any(|b| *b) {
        return false;
    }
    let valid: Vec<bool> = invalid.iter().map(|b| !b).collect();
    let (vidx, _) = edt_with_index(&valid, w, h);
    for (i, bad) in invalid.iter().enumerate() {
        if *bad {
            dem[i] = dem[vidx[i] as usize];
        }
    }
    true
}

/// 双线性上采样(同投影, 源网格从原点对齐)
fn upsample(src: &[f32], sw: usize, sh: usize, src_res: f64, dst_w: usize, dst_h: usize, dst_res: f64) -> Vec<f32> {
    let mut dst = vec![0f32; dst_w * dst_h];
    dst.par_chunks_mut(dst_w).enumerate().for_each(|(y, row)| {
        let gy = (y as f64 + 0.5) * dst_res / src_res - 0.5;
        let y0f = gy.floor();
        let fy = (gy - y0f) as f32;
        // 分辨率非整除时目标网格可略超源覆盖, 采样索引必须双向钳到边缘
        let y0 = (y0f.max(0.0) as usize).min(sh - 1);
        let y1 = (y0 + 1).min(sh - 1);
        for (x, out) in row.iter_mut().enumerate() {
            let gx = (x as f64 + 0.5) * dst_res / src_res - 0.5;
            let x0f = gx.floor();
            let fx = (gx - x0f) as f32;
            let x0 = (x0f.max(0.0) as usize).min(sw - 1);
            let x1 = (x0 + 1).min(sw - 1);
            let v00 = src[y0 * sw + x0];
            let v01 = src[y0 * sw + x1];
            let v10 = src[y1 * sw + x0];
            let v11 = src[y1 * sw + x1];
            *out = v00 * (1.0 - fx) * (1.0 - fy)
                + v01 * fx * (1.0 - fy)
                + v10 * (1.0 - fx) * fy
                + v11 * fx * fy;
        }
    });
    dst
}

fn levels_lut() -> [u8; 6] {
    [6, 5, 4, 3, 2, 1]
}

/// 升序分位数
pub fn pct_of(sorted: &[f32], p: f64) -> f32 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

/// 填充坝子内部完全包围的非坝子孔洞(4 邻域连通到边界的保持背景)
pub fn fill_interior_holes(basin: &mut [bool], w: usize, h: usize) {
    let n = w * h;
    let inv: Vec<bool> = basin.iter().map(|b| !b).collect();
    let mut lab = vec![0i32; n];
    let mut cur = 0i32;
    let mut touches = Vec::new();
    let mut stack: Vec<usize> = Vec::with_capacity(1024);
    for s0 in 0..n {
        if !inv[s0] || lab[s0] != 0 {
            continue;
        }
        cur += 1;
        lab[s0] = cur;
        stack.push(s0);
        let mut border = false;
        while let Some(i) = stack.pop() {
            let x = i % w;
            if x == 0 || x == w - 1 || i < w || i >= n - w {
                border = true;
            }
            let (cx, cy) = (x, i / w);
            for (dx, dy) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                let nx = cx as isize + dx;
                let ny = cy as isize + dy;
                if nx < 0 || ny < 0 || nx >= w as isize || ny >= h as isize {
                    continue;
                }
                let j = ny as usize * w + nx as usize;
                if inv[j] && lab[j] == 0 {
                    lab[j] = cur;
                    stack.push(j);
                }
            }
        }
        touches.push(border);
    }
    for i in 0..n {
        let l = lab[i];
        if l > 0 && !touches[(l - 1) as usize] {
            basin[i] = true;
        }
    }
}

/// 腐蚀(EDT 圆盘): 保留到背景距离 >= r 的像元
pub fn erode_round(mask: &mut [bool], w: usize, h: usize, r_px: f64) {
    if r_px < 1.0 {
        return;
    }
    let inv: Vec<bool> = mask.iter().map(|b| !b).collect();
    let (_, d) = edt_with_index(&inv, w, h);
    for (m, &v) in mask.iter_mut().zip(d.iter()) {
        *m = v >= r_px as f32;
    }
}

/// 膨胀(EDT 圆盘): 到集合(mask=true)距离 <= r 的像元置真
pub fn dilate_round(mask: &mut [bool], w: usize, h: usize, r_px: f64) {
    if r_px < 1.0 {
        return;
    }
    // EDT 源必须是集合本身: 返回每像元到最近集合像元的距离
    let (_, d) = edt_with_index(mask, w, h);
    for (m, &v) in mask.iter_mut().zip(d.iter()) {
        *m = v <= r_px as f32;
    }
}

/// 闭运算(先膨胀填缺口, 再腐蚀恢复边界)
pub fn closing_round(mask: &mut [bool], w: usize, h: usize, r_m: f64, res: f64) {
    let r = (r_m / res).round();
    if r < 1.0 {
        return;
    }
    dilate_round(mask, w, h, r);
    erode_round(mask, w, h, r);
}

/// 全流程运行。`cancel`: 返回 true 时中止。
pub fn run(
    params: &Params,
    progress: &dyn Fn(Progress) -> bool,
    cancelled: &std::sync::atomic::AtomicBool,
) -> Result<Outputs> {
    use std::sync::atomic::Ordering;
    let say = |stage: &str, pct: f32, msg: &str| -> Result<()> {
        if cancelled.load(Ordering::Relaxed) {
            return Err(CoreError::Cancelled);
        }
        progress(Progress { stage: stage.into(), pct, msg: msg.into() });
        Ok(())
    };

    // ---------- 1. 载入 5m DEM ----------
    say("载入", 0.0, "读取 DEM...")?;
    let (mut dem5, meta5) = geotiff::read_f32(&params.dem_path)?;
    let (w5, h5) = (meta5.width as usize, meta5.height as usize);
    let res_f = meta5.resolution();
    fill_nodata(&mut dem5, w5, h5);
    say("载入", 2.0, &format!("DEM {}x{} @{}m", w5, h5, res_f))?;

    // ---------- 2. 粗层 ----------
    say("粗层", 4.0, "生成中间分析层...")?;
    let (mut dem_c, cw, ch) = box_downsample(&dem5, w5, h5, res_f, params.coarse_res);
    let res_c = params.coarse_res;
    fill_nodata(&mut dem_c, cw, ch);

    // ---------- 3. 全填洼(坡位 TPI 专用; 坝子判据用原始地形) ----------
    say("填洼", 8.0, "Priority-Flood 全填洼(坡位用)...")?;
    let fr = fill_and_route(&dem_c, cw, ch, 99999.0);
    let filled_c = fr.filled.clone();

    // ---------- 4. 粗层因子(亚类起伏度 + 六级坡位) ----------
    say("因子", 18.0, "亚类起伏度 + geomorphon 六级坡位...")?;
    let relief2k_c = focal_relief(&filled_c, cw, ch, (params.relief_subclass_win / res_c) as usize);

    // 六级坡位(geomorphon, 对齐 GRASS r.geomorphon; 5m 原生分辨率计算):
    // 线视原理要求原始分辨率——窄谷(50-100m 宽)在粗层只有 2-4 像元,
    // 线视跨越后谷/坡下信号完全丢失(实测 25m 层谷/坡下为 0)。
    // rayon 并行; P 值 = 更高方向数(0..8), 越大越凹越近谷。
    say("坡位", 30.0, "geomorphon 5m 原生计算...")?;
    let search_px = params.slope_search_m / res_f;
    let skip_px = params.slope_skip_m / res_f;
    let flat_deg = params.slope_flat_deg;
    let pvals: Vec<u8> = (0..h5)
        .into_par_iter()
        .flat_map_iter(|y| {
            let mut row = vec![0u8; w5];
            if y > 0 && y < h5 - 1 {
                for (x, cell) in row.iter_mut().enumerate().take(w5 - 1).skip(1) {
                    let pat = crate::geomorphon::geomorphon_pattern(
                        &dem5, w5, h5, res_f, x, y, search_px, skip_px, flat_deg,
                    );
                    *cell = pat.count_higher() as u8;
                }
            }
            row
        })
        .collect::<Vec<u8>>();
    // ---------- 6. 上采样至成品分辨率 ----------
    say("上采样", 46.0, "坡位与起伏度插值...")?;
    let relief2k5 = upsample(&relief2k_c, cw, ch, res_c, w5, h5, res_f);
    drop(relief2k_c);
    // P 值直方图匹配到 legacy 六级分布(山脊 8.9 -> 坡上 14.4 -> 坡中 44.6
    // -> 平坡 3.5 -> 坡下 20.8 -> 山谷 7.8; P 越大越凹越近谷):
    let mut slope_pos = vec![0u8; w5 * h5];
    {
        let mut hist = [0u64; 9];
        for &v in pvals.iter() {
            hist[(v as usize).min(8)] += 1;
        }
        let total = (w5 * h5) as f64;
        let targets = [8.9f64, 14.4, 44.6, 3.5, 20.8, 7.8];
        let mut cuts = [8usize; 5]; // 各级上界(P 值), 缺省 8(兜底)
        let mut cum = 0f64;
        let mut ci = 0usize;
        for (pv, &h) in hist.iter().enumerate() {
            cum += h as f64;
            while ci < 5 && (cum / total) * 100.0 >= targets[..=ci].iter().sum::<f64>() {
                cuts[ci] = pv;
                ci += 1;
            }
        }
        for (out, &pv) in slope_pos.iter_mut().zip(pvals.iter()) {
            let p = (pv as usize).min(8);
            let mut lv = 1u8;
            for (ci, &cut) in cuts.iter().enumerate() {
                if p <= cut {
                    lv = levels_lut()[ci];
                    break;
                }
            }
            *out = lv;
        }
        say("坡位", 44.0, &format!("P 值切分点 {:?} (直方图匹配 legacy 分布)", cuts))?;
    }

    // 坡位后处理(对齐 notebook: 众数滤波 8 邻域 + 小斑(<200 像元@5m)蚕食)
    {
        let classes: [u8; 6] = [1, 2, 3, 4, 5, 6];
        for _ in 0..1 {
            mode_filter_pass(&mut slope_pos, w5, h5, &classes);
        }
        let min_cells = (params.slope_min_patch_m2 / (res_f * res_f)).ceil() as u64;
        let mut small = vec![false; w5 * h5];
        for c in classes {
            let m: Vec<bool> = slope_pos.iter().map(|&v| v == c).collect();
            let (lab, ln) = label_connected(&m, w5, h5);
            let mut sizes = vec![0u64; ln + 1];
            for &l in lab.iter() {
                if l > 0 {
                    sizes[l as usize] += 1;
                }
            }
            for i in 0..w5 * h5 {
                let l = lab[i];
                if l > 0 && sizes[l as usize] < min_cells {
                    small[i] = true;
                }
            }
        }
        if small.iter().any(|b| *b) {
            let src_ok: Vec<bool> = small.iter().map(|b| !b).collect();
            let (sidx, _) = edt_with_index(&src_ok, w5, h5);
            for (i, sm) in small.iter().enumerate() {
                if *sm {
                    slope_pos[i] = slope_pos[sidx[i] as usize];
                }
            }
        }
    }

    // ---------- 7/8. 坝子(OBIA 对象级判据; 原始地形粗层) ----------
    // 平坦判据全部基于原始 DEM: 无限制填洼的假湖面(深填区)保持真实
    // 起伏, 不会再伪造平坦骗过判据(Manba 定位的关键缺陷)。
    say("坝子", 58.0, "坝子识别(OBIA 对象级, 原始地形)...")?;
    drop(filled_c);
    let slope_c = slope_horn_degrees(&dem_c, cw, ch, res_c);
    let relief_c = focal_relief(&dem_c, cw, ch, 5); // 5x5 起伏
    // 区域高程百分位(2km 窗, 25m 粗层)
    let pctl_win = (params.basin_pctl_win_m / res_c) as usize;
    let mean_pctl = focal_mean(&dem_c, cw, ch, pctl_win);
    let th_pctl = params.basin_pctl_max; // 30 = 邻域内低于 70% 像元
    let th_slope = params.basin_slope_th; // 6°
    let th_relief = params.basin_relief_m; // 5m
    let n_c = cw * ch;
    let mut basin_c = vec![false; n_c];
    let mut n_cand = 0u64;
    for i in 0..n_c {

        // 简化百分位近似: 用相对偏差替代精确百分位(性能考虑)
        let is_low = dem_c[i] < mean_pctl[i] * (1.0 - th_pctl as f32 / 100.0 * 0.06); // 邻域均值偏移 2% 即视为低地
        if slope_c[i] < th_slope as f32 && relief_c[i] < th_relief as f32 && is_low {
            basin_c[i] = true;
            n_cand += 1;
        }
    }
    drop(mean_pctl);
    drop(slope_c);
    drop(relief_c);
    say("坝子", 59.0, &format!("平坦+低洼候选 {:.1}% ({} 像元)",
        100.0 * n_cand as f64 / n_c as f64, n_cand))?;
    // 上采样 5m
    let mut basin = vec![false; w5 * h5];
    let scale = (res_c / res_f).round() as usize;
    for y in 0..h5 {
        let sy = (y / scale).min(ch - 1);
        for x in 0..w5 {
            basin[y * w5 + x] = basin_c[sy * cw + (x / scale).min(cw - 1)];
        }
    }
    drop(basin_c);
    // 开运算(剔窄谷丝, 30m)
    let r_open = (30.0 / res_f).round();
    erode_round(&mut basin, w5, h5, r_open);
    dilate_round(&mut basin, w5, h5, r_open);
    // 面积过滤
    let (blab, bn) = label_connected(&basin, w5, h5);
    let mut bsizes = vec![0u64; bn + 1];
    for &l in blab.iter() {
        if l > 0 {
            bsizes[l as usize] += 1;
        }
    }
    let bmin = (params.basin_min_area_m2 / (res_f * res_f)) as u64;
    basin = (0..w5 * h5)
        .map(|i| blab[i] > 0 && bsizes[blab[i] as usize] >= bmin)
        .collect();
    // 填洞
    fill_interior_holes(&mut basin, w5, h5);
    // 平滑
    let r_sm = params.basin_smooth_m / res_f;
    dilate_round(&mut basin, w5, h5, r_sm);
    erode_round(&mut basin, w5, h5, r_sm);
    let n_basin = basin.iter().filter(|b| **b).count();
    say("坝子", 62.0, &format!("坝子 {:.2}% ({} 像元)",
        100.0 * n_basin as f64 / w5 as f64 / h5 as f64, n_basin))?;
    // ---------- 9. 精细亚类(5m 产品) ----------
    let mut sub5 = vec![0u8; w5 * h5];
    for i in 0..w5 * h5 {
        let z = dem5[i];
        sub5[i] = if (z as f64) < params.hill_z_max {
            if (relief2k5[i] as f64) < params.relief_low_hill { 1 } else { 2 }
        } else if z < 1000.0f32 {
            3
        } else if z < 3500.0f32 {
            4
        } else if z < 5000.0f32 {
            5
        } else {
            6
        };
    }
    for (i, b) in basin.iter().enumerate() {
        if *b {
            sub5[i] = 7;
        }
    }
    drop(relief2k5);

    // ---------- 10. 三图叠加: 六级坡位 x 坝子 x 海拔分区(calc_dxbw 规则) ----------
    say("坡位", 68.0, "六级坡位 + 三图叠加...")?;
    // E 四级海拔分区: 1<500m 2:500-800 3:800-1200 4:>=1200
    let mut zone = vec![0u8; w5 * h5];
    for i in 0..w5 * h5 {
        let z = dem5[i];
        zone[i] = if (z as f64) < params.hill_z_max {
            1
        } else if z < 800.0 {
            2
        } else if z < 1200.0 {
            3
        } else {
            4
        };
    }
    // 叠加(calc_dxbw 精确规则): 坝子优先, 坝内坡中上降级坡下/丘陵下部
    let mut terrain = vec![0u8; w5 * h5];
    for i in 0..w5 * h5 {
        let s = slope_pos[i];
        let e = zone[i];
        if basin[i] {
            if s <= 3 {
                terrain[i] = 1; // 山间/宽谷盆地
            } else if e == 1 {
                terrain[i] = 5; // 丘陵下部
            } else {
                terrain[i] = 8; // 山地坡下
            }
        } else if e == 1 {
            // 丘陵(<500m): 坡上=S5,6 坡中=S3,4 坡下=S1,2
            terrain[i] = match s {
                5..=6 => 3,
                3..=4 => 4,
                _ => 5,
            };
        } else {
            // 山地(>=500m): 坡上=S5,6 坡中=S3,4 坡下=S1,2
            terrain[i] = match s {
                5..=6 => 6,
                3..=4 => 7,
                _ => 8,
            };
        }
    }
    drop(dem5);
    // basin 留存至输出段(basin_mask.tif)
    // ---------- 11. 后处理(众数滤波 + 小图斑; 坝子保护) ----------
    // 与基准一致: 坝子(1)不在处理类别中 -> 形态只由核心化决定,
    // 不参与众数投票/不被去斑/不作为填充源, 杜绝窄脖颈经后处理回渗成坝子。
    say("后处理", 74.0, "众数滤波 + 小图斑去除...")?;
    let classes: [u8; 6] = [3, 4, 5, 6, 7, 8];
    for _ in 0..params.mode_filter_iter {
        mode_filter_pass(&mut terrain, w5, h5, &classes);
    }
    // 小图斑去除
    for _pass in 0..1 {
        let mut small = vec![false; w5 * h5];
        for &c in &classes {
            let m: Vec<bool> = terrain.iter().map(|&v| v == c).collect();
            let (lab, n) = label_connected(&m, w5, h5);
            let mut sizes = vec![0u64; n + 1];
            for &l in lab.iter() {
                if l > 0 {
                    sizes[l as usize] += 1;
                }
            }
            let th = (params.min_patch_m2 / (res_f * res_f)).ceil() as u64;
            for i in 0..w5 * h5 {
                let l = lab[i];
                if l > 0 && sizes[l as usize] < th {
                    small[i] = true;
                }
            }
        }
        if small.iter().any(|b| *b) {
            // 最近邻类别填充(EDT); 填充源排除坝子 -> 小斑只从坡位类取值
            let src_ok: Vec<bool> = small
                .iter()
                .zip(terrain.iter())
                .map(|(&s, &c)| !s && c != 1)
                .collect();
            let (sidx, _) = edt_with_index(&src_ok, w5, h5);
            for (i, s) in small.iter().enumerate() {
                if *s {
                    terrain[i] = terrain[sidx[i] as usize];
                }
            }
        }
    }

    // 盆地内小碎斑吞并: 山间盆地(编码1)内部、与外部不连通且面积
    // < 10000m² 的非盆地碎斑(细碎坡下/坡中)整体并入盆地
    // ——盆地应为整体连片的平坦区(Manba 口径, 单斑上限 10000m²)
    {
        let hole_max = (10_000.0 / (res_f * res_f)).ceil() as u64; // 400 像元@5m
        let holes: Vec<bool> = terrain
            .iter()
            .map(|&t| t != 1 && t != 0)
            .collect();
        // 4 邻域连通的"非盆地"区块; 与图像边界接触=通向外部(保留)
        let mut lab = vec![0i32; w5 * h5];
        let mut cur = 0i32;
        let mut sizes: Vec<u64> = Vec::new();
        let mut touches = Vec::new();
        let mut stack: Vec<usize> = Vec::with_capacity(1024);
        for s0 in 0..w5 * h5 {
            if !holes[s0] || lab[s0] != 0 {
                continue;
            }
            cur += 1;
            lab[s0] = cur;
            stack.push(s0);
            let mut sz = 0u64;
            let mut border = false;
            while let Some(i) = stack.pop() {
                sz += 1;
                let x = i % w5;
                if x == 0 || x == w5 - 1 || i < w5 || i >= w5 * h5 - w5 {
                    border = true;
                }
                let (cx, cy) = (x, i / w5);
                for (dx, dy) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                    let nx = cx as isize + dx;
                    let ny = cy as isize + dy;
                    if nx < 0 || ny < 0 || nx >= w5 as isize || ny >= h5 as isize {
                        continue;
                    }
                    let j = ny as usize * w5 + nx as usize;
                    if holes[j] && lab[j] == 0 {
                        lab[j] = cur;
                        stack.push(j);
                    }
                }
            }
            sizes.push(sz);
            touches.push(border);
        }
        let mut merged = 0u64;
        for i in 0..w5 * h5 {
            let l = lab[i];
            if l > 0 && !touches[(l - 1) as usize] && sizes[(l - 1) as usize] <= hole_max {
                terrain[i] = 1;
                merged += 1;
            }
        }
        drop(lab);
        say("盆地", 70.0, &format!("盆地内碎斑吞并 {} 像元 (单斑<={}m²)",
            merged, 10_000))?;
    }

    // ---------- 12. 输出 ----------
    say("输出", 88.0, "写栅格与报告...")?;
    let out_dir = Path::new(&params.out_dir);
    std::fs::create_dir_all(out_dir)?;
    let terrain_path = out_dir.join("terrain_position.tif");
    let mut cmap = [[0u8; 3]; 256];
    cmap[1] = [51, 178, 229];
    cmap[3] = [250, 217, 89];
    cmap[4] = [217, 237, 166];
    cmap[5] = [153, 199, 102];
    cmap[6] = [250, 165, 60];
    cmap[7] = [222, 100, 50];
    cmap[8] = [107, 68, 35];
    geotiff::write_u8_cmap(&terrain_path, &meta5, &terrain, &cmap)?;
    let sub_path = out_dir.join("geomorph_subclass.tif");
    let mut sub_cmap = [[0u8; 3]; 256];
    sub_cmap[1] = [180, 230, 150];
    sub_cmap[2] = [110, 195, 110];
    sub_cmap[3] = [250, 225, 130];
    sub_cmap[4] = [235, 170, 90];
    sub_cmap[5] = [205, 110, 75];
    sub_cmap[6] = [150, 65, 60];
    sub_cmap[7] = [85, 185, 235];
    geotiff::write_u8_cmap(&sub_path, &meta5, &sub5, &sub_cmap)?;
    // 三张叠加中间成果: 坡位图 / 坝子图 / 丘陵山地图(可独立核查)
    let mut pos_cmap = [[0u8; 3]; 256];
    pos_cmap[1] = [168, 112, 72];  // 坡上
    pos_cmap[2] = [222, 196, 120]; // 坡中
    pos_cmap[3] = [132, 168, 96];  // 坡下
    let _ = geotiff::write_u8_cmap(out_dir.join("slope_position.tif"), &meta5, &slope_pos, &pos_cmap);
    drop(slope_pos);
    let mut basin_cmap = [[0u8; 3]; 256];
    basin_cmap[1] = [51, 178, 229]; // 坝子
    let basin_u8: Vec<u8> = basin.iter().map(|&b| b as u8).collect();
    let _ = geotiff::write_u8_cmap(out_dir.join("basin_mask.tif"), &meta5, &basin_u8, &basin_cmap);
    drop(basin_u8);
    drop(basin);
    let mut zone_cmap = [[0u8; 3]; 256];
    zone_cmap[1] = [180, 230, 150]; // 丘陵
    zone_cmap[2] = [235, 170, 90];  // 山地
    let _ = geotiff::write_u8_cmap(out_dir.join("hill_mountain_zone.tif"), &meta5, &zone, &zone_cmap);
    drop(zone);


    let names = [
        (1u8, "山间盆地"),
        (3, "丘陵上部"),
        (4, "丘陵中部"),
        (5, "丘陵下部"),
        (6, "山地坡上"),
        (7, "山地坡中"),
        (8, "山地坡下"),
    ];
    // 全分辨率统计(主 DEM 像元 25m²)
    let stat_order: [(u8, &str); 7] = [
        (1, "山间盆地"),
        (6, "山地坡上"),
        (7, "山地坡中"),
        (8, "山地坡下"),
        (3, "丘陵上部"),
        (4, "丘陵中部"),
        (5, "丘陵下部"),
    ];
    let stats: Vec<(u8, f64)> = stat_order
        .iter()
        .map(|&(c, _)| {
            (c, terrain.iter().filter(|&&v| v == c).count() as f64 * res_f * res_f / 1e6)
        })
        .collect();
    let mut lines = vec!["地形部位划分面积统计".to_string(), "=".repeat(46)];
    let mut total = 0f64;
    for (v, name) in names {
        let c = terrain.iter().filter(|&&x| x == v).count() as f64;
        total += c;
        lines.push(format!("{}  {:<8} {:>10.2} km²  {:>6.2}%", v, name, c * res_f * res_f / 1e6, 100.0 * c / (terrain.len() as f64)));
    }
    let _ = total;
    let report = lines.join("\n");
    std::fs::write(out_dir.join("class_report.txt"), &report)?;

    Ok(Outputs {
        terrain,
        subclass: sub5,
        meta5,
        report,
        stats,
    })
}

fn mode_filter_pass(arr: &mut [u8], w: usize, h: usize, classes: &[u8]) {
    let src = arr.to_vec();
    let in_domain: Vec<bool> = src.iter().map(|v| classes.contains(v)).collect();
    let out: Vec<u8> = (0..w * h).into_par_iter().map(|i| {
        // 保护类(如坝子)原值保持
        if !in_domain[i] {
            return src[i];
        }
        let x = i % w;
        let y = i / w;
        let mut cnt = [0u16; 16];
        for dy in -1isize..=1 {
            for dx in -1isize..=1 {
                let yy = (y as isize + dy).clamp(0, h as isize - 1) as usize;
                let xx = (x as isize + dx).clamp(0, w as isize - 1) as usize;
                let j = yy * w + xx;
                // 投票只计参与处理的类别, 保护类不计票
                if in_domain[j] {
                    cnt[src[j] as usize] += 1;
                }
            }
        }
        // 众数(同票保留原值)
        let mut best = src[i];
        let mut best_n = cnt[src[i] as usize];
        for &c in classes {
            if cnt[c as usize] > best_n {
                best_n = cnt[c as usize];
                best = c;
            }
        }
        best
    }).collect();
    arr.copy_from_slice(&out);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dilate_direction() {
        // 回归: dilate 曾误用"到背景距离"导致全图饱和
        let (w, h) = (50usize, 50usize);
        let mut m = vec![false; w * h];
        m[25 * w + 25] = true;
        dilate_round(&mut m, w, h, 3.0);
        assert!(m[25 * w + 25], "源保持");
        assert!(m[25 * w + 28], "右侧 3px 内膨胀");
        assert!(!m[28 * w + 28], "对角 sqrt(18)>3 不膨胀");
        // 远角不受影响
        assert!(!m[0]);
    }

    #[test]
    fn close_round_monotonic() {
        // 闭运算单调性: X ⊆ closing(X), 像元数不得减少
        let (w, h) = (800usize, 600usize);
        let mut m = vec![false; w * h];
        // 伪随机带状图案 ~35% 覆盖
        let mut seed = 12345u64;
        for i in 0..w * h {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            m[i] = (seed >> 33) % 100 < 35;
        }
        let before = m.iter().filter(|b| **b).count();
        closing_round(&mut m, w, h, 4.0, 1.0); // r=4px
        let after = m.iter().filter(|b| **b).count();
        assert!(after >= before, "闭运算收缩: before={before} after={after}");
    }

    #[test]
    fn upsample_non_integral_ratio_no_oob() {
        // 复现线上闪退: 5m→25m 非整除(15813*0.2=3162.6, 粗层取整 3162),
        // 目标网格最后一列/行的采样索引曾越界 1 个像元
        let (cw, ch) = (3162usize, 1900usize);
        let src = vec![1.0f32; cw * ch];
        let dst = upsample(&src, cw, ch, 25.0, 15813, 9500, 5.0);
        assert_eq!(dst.len(), 15813 * 9500);
        assert!(dst.iter().all(|&v| (v - 1.0).abs() < 1e-6));
    }

    #[test]
    fn upsample_identity() {
        let src = vec![2.5f32; 25];
        let dst = upsample(&src, 5, 5, 25.0, 5, 5, 25.0);
        assert!(dst.iter().all(|&v| v == 2.5));
    }
}
