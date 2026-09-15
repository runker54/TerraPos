//! 盆地对象识别: 低平候选 -> 宽度核心 -> 对象检验 -> 受约束边界重建。
//!
//! 候选判据全部使用米制阈值并在所属模块边界一次换算; 对象接受采用
//! 面积/宽度/内起伏/围限/水文连通/多尺度保持多重判据; 通过对象仅从
//! 宽度核心在原候选域内重建完整边界, 不输出纯腐蚀核心。

use crate::distance::edt_with_index;
use crate::error::Result;
use crate::geomorphon::MorphEvidence;
use crate::hydro::HydroModel;
use crate::input::RasterShape;
use crate::pipeline::BasinTendency;
use crate::scale::ScalePyramid;
use crate::slope_unit::{SlopeGeometry, SlopeUnits};

/// 盆地识别配置
#[derive(Debug, Clone)]
pub struct BasinConfig {
    pub tendency: BasinTendency,
    pub min_area_m2: f64,
}

/// 单个盆地对象诊断特征
#[derive(Debug, Clone)]
pub struct BasinObject {
    pub id: u32,
    pub area_m2: f64,
    pub max_width_m: f32,
    pub median_width_m: f32,
    pub inner_relief_m: f32,
    pub surround_rise_m: f32,
    pub stream_connected: bool,
    pub closed_depression: bool,
    pub accepted: bool,
}

/// 盆地识别产品: 候选/核心/重建掩膜与对象表
#[derive(Debug, Clone, Default)]
pub struct BasinResult {
    pub candidate: Vec<bool>,
    pub core: Vec<bool>,
    pub mask: Vec<bool>,
    pub objects: Vec<BasinObject>,
}

/// 亚类坡度上限(低丘/高丘/低山/中山/高山/极高山)
const SLOPE_LIMITS: [f32; 6] = [6.0, 5.5, 5.0, 5.0, 4.0, 4.0];

/// 亚类判定(与 scale::scale_bounds_for 同口径: <500m 丘陵按起伏 200 分低/高丘)
fn subclass_index(elev: f32, relief_m: f32) -> usize {
    if elev < 500.0 {
        if relief_m < 200.0 {
            0
        } else {
            1
        }
    } else if elev < 1000.0 {
        2
    } else if elev < 3500.0 {
        3
    } else if elev < 5000.0 {
        4
    } else {
        5
    }
}

/// 组件内参考坡面宽度 L 的中位数(米)
fn median_l_for(comp: &[usize], geometry: &SlopeGeometry) -> f32 {
    let mut ls: Vec<f32> = comp
        .iter()
        .map(|&i| geometry.local_width_m[i])
        .filter(|v| v.is_finite())
        .collect();
    ls.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if ls.is_empty() {
        500.0
    } else {
        ls[ls.len() / 2]
    }
}

/// 在组件范围内填充掩膜内部小洞(接触组件外部的保持不填)
fn fill_small_holes(mask: &mut [bool], w: usize, h: usize, comp: &[usize], hole_max_cells: usize) {
    let n = w * h;
    let mut lab = vec![0u32; n];
    let mut cur = 0u32;
    let mut sizes: Vec<usize> = Vec::new();
    let mut touches_out: Vec<bool> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    for &s0 in comp {
        if mask[s0] || lab[s0] != 0 {
            continue;
        }
        cur += 1;
        lab[s0] = cur;
        stack.push(s0);
        let mut sz = 0usize;
        let mut touches = false;
        while let Some(i) = stack.pop() {
            sz += 1;
            let x = (i % w) as i64;
            let y = (i / w) as i64;
            for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    touches = true;
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if !mask[j] && lab[j] == 0 {
                    lab[j] = cur;
                    stack.push(j);
                }
            }
        }
        sizes.push(sz);
        touches_out.push(touches);
    }
    for i in 0..n {
        let l = lab[i];
        if l > 0 && !touches_out[(l - 1) as usize] && sizes[(l - 1) as usize] <= hole_max_cells {
            mask[i] = true;
        }
    }
}

/// 检测盆地对象并重建完整边界
#[allow(clippy::too_many_arguments)]
pub fn detect_basins(
    dem: &[f32],
    valid: &[bool],
    shape: RasterShape,
    hydro: &HydroModel,
    units: &SlopeUnits,
    geometry: &SlopeGeometry,
    pyramid: &ScalePyramid,
    morph: &MorphEvidence,
    cfg: &BasinConfig,
) -> Result<BasinResult> {
    let w = shape.width;
    let h = shape.height;
    let n = w * h;
    let res = shape.resolution_m;
    let cell_area = res * res;
    let tendency_factor = match cfg.tendency {
        BasinTendency::Strict => 0.85f32,
        BasinTendency::Standard => 1.00,
        BasinTendency::Loose => 1.15,
    };

    // 闭合洼地汇流区(规格 8.3: 独立水文子系统): 从深洼环沿 flow_to
    // 反向 BFS; 洼地豁免只授予区内像元, 深洼环(溢口边缘)本身不豁免
    let mut closed_region = vec![false; n];
    {
        let mut stack: Vec<usize> = (0..n).filter(|&i| hydro.deep_sink[i]).collect();
        while let Some(i) = stack.pop() {
            let x = (i % w) as i64;
            let y = (i / w) as i64;
            for (dx, dy) in [
                (-1i64, 0i64),
                (1, 0),
                (0, -1),
                (0, 1),
                (-1, -1),
                (1, -1),
                (-1, 1),
                (1, 1),
            ] {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if !closed_region[j] && hydro.flow_to[j] == i as u32 {
                    closed_region[j] = true;
                    stack.push(j);
                }
            }
        }
    }

    // 1) 像元候选: 下部证据 + 低 q + 低坡 + 微平 + 近水/洼地 + 多尺度保持
    let mut candidate = vec![false; n];
    let mut relief_micro: Vec<f32> = Vec::new();
    for i in 0..n {
        if !valid[i] || !morph.lower[i].is_finite() {
            continue;
        }
        if morph.lower[i] <= morph.upper[i] {
            continue;
        }
        let q = geometry.relative_position[i];
        if !q.is_finite() || q >= 0.45 {
            continue;
        }
        let h_star = pyramid.characteristic_relief_m[i].max(0.0);
        let l = geometry.local_width_m[i];
        let l_ref = if l.is_finite() {
            l
        } else {
            2.0 * pyramid.characteristic_scale_m[i]
        };
        let slope_limit = SLOPE_LIMITS[subclass_index(dem[i], h_star)] * tendency_factor;
        if morph.slope_deg[i] > slope_limit {
            continue;
        }
        let h_micro_max = (0.03 * h_star).clamp(3.0, 12.0) * tendency_factor;
        if relief_micro.is_empty() {
            let r_micro_m = (0.05 * l_ref.max(500.0)).clamp(25.0, 100.0);
            let r_micro_px = ((r_micro_m as f64 / res).round() as usize).max(1);
            relief_micro = crate::filter::focal_relief(dem, w, h, 2 * r_micro_px + 1);
        }
        if relief_micro[i] > h_micro_max {
            continue;
        }
        let hand_max = (0.02 * h_star + 2.0).clamp(5.0, 20.0) * tendency_factor;
        let hand_i = geometry.hand_m[i];
        let in_stream_band = hand_i.is_finite() && hand_i <= hand_max;
        // 豁免仅授予无水文锚(链上无河网, HAND 为 NaN)的洼内像元
        // (喀斯特型洼底); 有锚的内坡像元仍按 HAND 判据执行
        let in_closed = closed_region[i] && !geometry.hand_m[i].is_finite();
        if !(in_stream_band || in_closed) {
            continue;
        }
        // 三连续尺度低平保持: 用逐尺度四分位距(对远处围限高地稳健)
        let mut persist = 0usize;
        let mut best = 0usize;
        for layer in &pyramid.layers {
            let iqr = layer.quartile_spread_m[i];
            if iqr.is_finite() && iqr <= h_micro_max.max(12.0) {
                persist += 1;
                best = best.max(persist);
            } else {
                persist = 0;
            }
        }
        if best < 3 && !(in_closed && best >= 1) {
            continue;
        }
        candidate[i] = true;
    }

    // 2) 对象组成: 同一谷底单元内的候选聚合为同一对象(规格 11.3),
    // 谷屏障上的候选(unit 0)并入相邻单元对象, 避免河链把盆底切碎
    let mut objects: Vec<BasinObject> = Vec::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    let outside: Vec<bool> = (0..n).map(|i| valid[i] && !candidate[i]).collect();
    let (_src, dist_edge) = edt_with_index(&outside, w, h);
    let mut lab = vec![usize::MAX; n];
    {
        let mut obj_of_unit: std::collections::HashMap<u32, usize> =
            std::collections::HashMap::new();
        for i in 0..n {
            if !candidate[i] || units.unit_id[i] == 0 {
                continue;
            }
            let u = units.unit_id[i];
            let oi = *obj_of_unit.entry(u).or_insert_with(|| {
                components.push(Vec::new());
                components.len() - 1
            });
            lab[i] = oi;
            components[oi].push(i);
        }
        // 屏障候选两遍赋值: 各自归入邻接单元对象, 屏障像元之间
        // 不得互为桥(否则河链会把地形上分隔的对象串联)
        let mut barrier_targets: Vec<(usize, Option<usize>)> = Vec::new();
        for i in 0..n {
            if !candidate[i] || units.unit_id[i] != 0 {
                continue;
            }
            let x = (i % w) as i64;
            let y = (i / w) as i64;
            let mut target: Option<usize> = None;
            for (dx, dy) in [
                (-1i64, 0i64),
                (1, 0),
                (0, -1),
                (0, 1),
                (-1, -1),
                (1, -1),
                (-1, 1),
                (1, 1),
            ] {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if candidate[j] && units.unit_id[j] != 0 && lab[j] != usize::MAX {
                    target = Some(match target {
                        Some(t) => t.min(lab[j]),
                        None => lab[j],
                    });
                }
            }
            barrier_targets.push((i, target));
        }
        for (i, target) in barrier_targets {
            match target {
                Some(oi) => {
                    lab[i] = oi;
                    components[oi].push(i);
                }
                None => {
                    lab[i] = components.len();
                    components.push(vec![i]);
                }
            }
        }
    }
    for (k, comp) in components.iter().enumerate() {
        if comp.is_empty() {
            continue;
        }
        let area_m2 = comp.len() as f64 * cell_area;
        // 宽度分布: 组件内到边界的 EDT
        let mut widths: Vec<f32> = Vec::with_capacity(comp.len());
        let mut max_width = 0f32;
        for &i in comp {
            let w_i = 2.0 * dist_edge[i] * res as f32;
            widths.push(w_i);
            if w_i > max_width {
                max_width = w_i;
            }
        }
        widths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median_width = widths[widths.len() / 2];
        let med_l = median_l_for(comp, geometry);
        let inner_med = {
            let mut zs: Vec<f32> = comp.iter().map(|&i| dem[i]).collect();
            zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            zs[zs.len() / 2]
        };
        let mut hs: Vec<f32> = comp
            .iter()
            .map(|&i| pyramid.characteristic_relief_m[i])
            .filter(|v| v.is_finite())
            .collect();
        hs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med_h = if hs.is_empty() { 50.0 } else { hs[hs.len() / 2] };
        let required_surround = 10.0f32.max(0.05 * med_h);
        // 方向围限剖面: 16 方向从等效半径向外扫描 ring 半径,
        // 统计剖面最大高差达标的方向占比(环带统计会被候选域
        // 形状稀释, 方向剖面才表达"周围存在围限高地"的语义)
        let ring_radius = (0.15 * med_l).clamp(150.0, 500.0);
        let cx_f = comp.iter().map(|&i| (i % w) as f64).sum::<f64>() / comp.len() as f64;
        let cy_f = comp.iter().map(|&i| (i / w) as f64).sum::<f64>() / comp.len() as f64;
        let r_eq_m = area_m2.sqrt().max(res);
        let r_start_px = (r_eq_m / res).ceil() as i64;
        // 扫描距离取 2*ring: 保证围限高地在近邻搜索范围内
        let r_end_px = ((r_eq_m + 2.0 * ring_radius as f64) / res).ceil() as i64;
        let sample_z = |fx: f64, fy: f64| -> f32 {
            let x0 = fx.floor() as isize;
            let y0 = fy.floor() as isize;
            let g = |xx: isize, yy: isize| -> f32 {
                if xx < 0 || yy < 0 || xx >= w as isize || yy >= h as isize {
                    f32::NAN
                } else {
                    dem[(yy as usize) * w + xx as usize]
                }
            };
            let z00 = g(x0, y0);
            let z10 = g(x0 + 1, y0);
            let z01 = g(x0, y0 + 1);
            let z11 = g(x0 + 1, y0 + 1);
            let gx = fx - x0 as f64;
            let gy = fy - y0 as f64;
            let v = z00 as f64 * (1.0 - gx) * (1.0 - gy)
                + z10 as f64 * gx * (1.0 - gy)
                + z01 as f64 * (1.0 - gx) * gy
                + z11 as f64 * gx * gy;
            if v.is_finite() { v as f32 } else { f32::NAN }
        };
        let mut dir_rises: Vec<f32> = Vec::new();
        for d in 0..16 {
            let ang = d as f64 * std::f64::consts::FRAC_PI_8;
            let (dx, dy) = (ang.cos(), ang.sin());
            let mut max_z = f32::NEG_INFINITY;
            for s in r_start_px..=r_end_px {
                let z = sample_z(
                    cx_f + dx * s as f64,
                    cy_f + dy * s as f64,
                );
                if z.is_finite() && z > max_z {
                    max_z = z;
                }
            }
            if max_z.is_finite() {
                dir_rises.push(max_z - inner_med);
            }
        }
        dir_rises.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let hit_ratio = if dir_rises.is_empty() {
            0.0
        } else {
            dir_rises
                .iter()
                .filter(|&&r| r >= required_surround)
                .count() as f64
                / dir_rises.len() as f64
        };
        let surround_rise = dir_rises.last().copied().unwrap_or(0.0);
        let surround_ok = hit_ratio >= 0.6;
        let stream_connected = comp.iter().any(|&i| hydro.stream_level[i] > 0);
        let closed_depression = comp.iter().any(|&i| closed_region[i]);
        let id = k as u32 + 1;
        // 对象判据(米制阈值一次生成)
        let w_min = (0.12 * med_l).clamp(150.0, 500.0);
        let a_min = cfg
            .min_area_m2
            .max((0.015 * med_l as f64 * med_l as f64).clamp(50_000.0, 500_000.0));
        let h_inner_max = (0.08 * med_h).clamp(10.0, 30.0);
        let has_core = max_width >= w_min && median_width >= 0.5 * w_min;
        // 内部起伏在宽度核心区计算(核心代表盆底平坦部分,
        // 全候选域会混入坡脚过渡带高值)
        let inner_relief = {
            let mut zs: Vec<f32> = comp
                .iter()
                .filter(|&&i| dist_edge[i] * 2.0 >= w_min)
                .map(|&i| dem[i])
                .collect();
            if zs.is_empty() {
                0.0
            } else {
                zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
                zs[(zs.len() as f32 * 0.95) as usize % zs.len()]
                    - zs[(zs.len() as f32 * 0.05) as usize % zs.len()]
            }
        };
        let accepted = area_m2 >= a_min
            && has_core
            && inner_relief <= h_inner_max
            && surround_ok
            && (stream_connected || closed_depression);
        objects.push(BasinObject {
            id,
            area_m2,
            max_width_m: max_width,
            median_width_m: median_width,
            inner_relief_m: inner_relief,
            surround_rise_m: surround_rise,
            stream_connected,
            closed_depression,
            accepted,
        });
    }

    // 3) 核心与受约束重建: 接受对象从宽度核心在原候选组件内 geodesic 扩展
    let mut core = vec![false; n];
    let mut mask = vec![false; n];
    for (k, obj) in objects.iter().enumerate() {
        if !obj.accepted {
            continue;
        }
        let comp = &components[k];
        let w_min = (0.12 * median_l_for(comp, geometry)).clamp(150.0, 500.0);
        let mut obj_core = Vec::new();
        for &i in comp {
            if dist_edge[i] * 2.0 * res as f32 >= w_min {
                core[i] = true;
                obj_core.push(i);
                mask[i] = true;
            }
        }
        if obj_core.is_empty() {
            continue;
        }
        let mut stack = obj_core.clone();
        while let Some(i) = stack.pop() {
            let x = (i % w) as i64;
            let y = (i / w) as i64;
            for (dx, dy) in [
                (-1i64, 0i64),
                (1, 0),
                (0, -1),
                (0, 1),
                (-1, -1),
                (1, -1),
                (-1, 1),
                (1, 1),
            ] {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if candidate[j] && lab[j] == k && !mask[j] {
                    mask[j] = true;
                    stack.push(j);
                }
            }
        }
        let hole_max_cells = ((0.1 * obj.area_m2).min(20_000.0) / cell_area).ceil() as usize;
        fill_small_holes(&mut mask, w, h, comp, hole_max_cells);
    }

    Ok(BasinResult {
        candidate,
        core,
        mask,
        objects,
    })
}
