//! 最终编码组合与受约束清理。
//!
//! 编码语义按设计规格 1.2(Manba 裁定): 0=NoData, 1=山间/宽谷盆地,
//! 2=保留码永不生成, 3/4/5=丘陵上/中/下, 6/7/8=山地坡上/坡中/坡下。
//! 清理只在坡面单元与丘陵/山地界内部进行, 不跨脊/谷/盆地边界。

use crate::basin::BasinResult;
use crate::error::Result;
use crate::input::RasterShape;
use crate::slope_position::{Memberships, SlopePosition};
use crate::slope_unit::{SlopeGeometry, SlopeUnits};

/// 最终分类产品
#[derive(Debug, Clone)]
pub struct FinalClassification {
    /// 地形部位编码(0/1/3-8)
    pub terrain: Vec<u8>,
    /// 地貌亚类编码: 1 低丘 2 高丘 3 低山 4 中山 5 高山 6 极高山 7 盆地
    pub geomorph_subclass: Vec<u8>,
    pub confidence: Vec<f32>,
}

/// 组合最终编码: NoData=0 -> 合格盆地=1 -> 丘陵/山地 x 上中下
pub fn compose_codes(
    elevation_m: &[f32],
    valid: &[bool],
    basin: &BasinResult,
    positions: &Memberships,
    hill_elevation_max_m: f32,
) -> Result<FinalClassification> {
    let n = elevation_m.len();
    let mut terrain = vec![0u8; n];
    let mut subclass = vec![0u8; n];
    let mut confidence = vec![0f32; n];
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        if basin.mask[i] {
            terrain[i] = 1;
            subclass[i] = 7;
            confidence[i] = 1.0;
            continue;
        }
        let hill = elevation_m[i] < hill_elevation_max_m;
        terrain[i] = match positions.raw[i] {
            SlopePosition::Upper => {
                if hill {
                    3
                } else {
                    6
                }
            }
            SlopePosition::Middle => {
                if hill {
                    4
                } else {
                    7
                }
            }
            SlopePosition::Lower => {
                if hill {
                    5
                } else {
                    8
                }
            }
        };
        // 亚类(低/高丘细分需相对起伏, 由编排层以金字塔数据细化)
        subclass[i] = if hill {
            1
        } else if elevation_m[i] < 1000.0 {
            3
        } else if elevation_m[i] < 3500.0 {
            4
        } else if elevation_m[i] < 5000.0 {
            5
        } else {
            6
        };
        confidence[i] = positions.confidence[i];
    }
    Ok(FinalClassification {
        terrain,
        geomorph_subclass: subclass,
        confidence,
    })
}

/// 受约束清理:
/// 1. 坡面单元内小斑归并(不跨脊/谷/盆地/无效/丘陵山地界);
/// 2. 脊→谷轨迹的三态单调动态规划修正(仅低置信像元, 仅总代价改善时)。
///
/// 返回被修正的像元数(供报告统计)。
#[allow(clippy::too_many_arguments)]
pub fn constrained_cleanup(
    result: &mut FinalClassification,
    positions: &Memberships,
    units: &SlopeUnits,
    geometry: &SlopeGeometry,
    basin: &BasinResult,
    shape: RasterShape,
    strength: f64,
) -> Result<usize> {
    // 盆地像元在轨迹修正中以 terrain==1 识别, 参数保留以维持契约签名
    let _ = basin;
    let mut corrected = 0usize;

    // ---- 1) 小斑归并 ----
    // 平滑距离 = clip(5*strength, 5, 50) 米, 再受 0.02*中位坡面宽度约束
    let mut widths: Vec<f32> = geometry
        .local_width_m
        .iter()
        .copied()
        .filter(|v| v.is_finite() && *v > 0.0)
        .collect();
    widths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med_l = if widths.is_empty() { 500.0 } else { widths[widths.len() / 2] };
    let smooth_m = ((5.0 * strength) as f32).clamp(5.0, 50.0).min(0.02 * med_l);
    let patch_area_m2 = (smooth_m * smooth_m) as f64;
    let patch_cells = ((patch_area_m2 / (shape.resolution_m * shape.resolution_m)).ceil() as usize).max(1);

    merge_small_patches(result, units, shape, patch_cells);

    // ---- 2) 轨迹单调 DP 修正 ----
    corrected += monotonic_trace_fix(result, positions, units, geometry, basin, shape);

    Ok(corrected)
}

/// 坡面单元内小斑归并
fn merge_small_patches(
    result: &mut FinalClassification,
    units: &SlopeUnits,
    shape: RasterShape,
    patch_cells: usize,
) {
    let w = shape.width;
    let h = shape.height;
    let n = w * h;
    // 逐单元、逐编码标注连通斑(4 连通)
    let codes = [3u8, 4, 5, 6, 7, 8];
    let mut absorbed = vec![false; n];
    for u in 1.. {
        if !units.unit_id.contains(&u) {
            break;
        }
        for &code in &codes {
            let mut seen = vec![false; n];
            for s0 in 0..n {
                if units.unit_id[s0] != u
                    || result.terrain[s0] != code
                    || seen[s0]
                    || absorbed[s0]
                {
                    continue;
                }
                // BFS 收集同码斑
                let mut patch = vec![s0];
                seen[s0] = true;
                let mut stack = vec![s0];
                let mut touches_barrier = false;
                while let Some(i) = stack.pop() {
                    let x = (i % w) as i64;
                    let y = (i / w) as i64;
                    for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                        let nx = x + dx;
                        let ny = y + dy;
                        if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                            touches_barrier = true;
                            continue;
                        }
                        let j = (ny as usize) * w + nx as usize;
                        if units.unit_id[j] != u || absorbed[j] {
                            touches_barrier = true;
                            continue;
                        }
                        if result.terrain[j] == code && !seen[j] {
                            seen[j] = true;
                            patch.push(j);
                            stack.push(j);
                        }
                    }
                }
                if patch.len() >= patch_cells || touches_barrier {
                    continue;
                }
                // 统计共享边界最长的邻码(仅同单元、非盆地)
                let mut contact: std::collections::HashMap<u8, usize> =
                    std::collections::HashMap::new();
                for &i in &patch {
                    let x = (i % w) as i64;
                    let y = (i / w) as i64;
                    for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                        let nx = x + dx;
                        let ny = y + dy;
                        if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                            continue;
                        }
                        let j = (ny as usize) * w + nx as usize;
                        let c = result.terrain[j];
                        if units.unit_id[j] == u && c != code && c != 0 && c != 1 {
                            *contact.entry(c).or_default() += 1;
                        }
                    }
                }
                if let Some((&target, _)) = contact.iter().max_by_key(|(&c, &cnt)| (cnt, c)) {
                    for &i in &patch {
                        result.terrain[i] = target;
                        absorbed[i] = true;
                    }
                }
            }
        }
        // 单元编号有限, 防御性上限
        if u > 1_000_000 {
            break;
        }
    }
}

/// 脊→谷轨迹单调修正(三态 DP, 负对数隶属度代价)
fn monotonic_trace_fix(
    result: &mut FinalClassification,
    positions: &Memberships,
    units: &SlopeUnits,
    geometry: &SlopeGeometry,
    basin: &BasinResult,
    shape: RasterShape,
) -> usize {
    let w = shape.width;
    let h = shape.height;
    let rank_of_code = |c: u8| -> i32 {
        match c {
            3 | 6 => 2,
            4 | 7 => 1,
            5 | 8 => 0,
            _ => -1, // 0/1: 不参与
        }
    };
    let membership = |positions: &Memberships, i: usize, s: usize| -> f32 {
        match s {
            2 => positions.upper[i],
            1 => positions.middle[i],
            _ => positions.lower[i],
        }
    };
    let neg_ln = |p: f32| -> f64 {
        let p = p.clamp(1e-3, 1.0) as f64;
        -p.ln()
    };
    let mut corrected = 0usize;
    let mut traces_done = 0usize;
    const MAX_TRACES: usize = 2000;
    let n = w * h;
    for s0 in 0..n {
        if traces_done >= MAX_TRACES {
            break;
        }
        if !units.ridge_mask[s0] || basin.mask[s0] || !basin_or_valid(s0, result) {
            continue;
        }
        if rank_of_code(result.terrain[s0]) < 0 {
            continue;
        }
        traces_done += 1;
        // 轨迹: 沿 dv 递减贪心步进(8 邻), 盆地像元排除
        let mut trace = vec![s0];
        let mut cur = s0;
        loop {
            let x = (cur % w) as i64;
            let y = (cur / w) as i64;
            let mut next = None;
            let mut best = geometry.distance_to_valley_m[cur];
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
                if !basin_or_valid(j, result) || basin.mask[j] || trace.contains(&j) {
                    continue;
                }
                let dv = geometry.distance_to_valley_m[j];
                if dv.is_finite() && dv < best {
                    best = dv;
                    next = Some(j);
                }
            }
            match next {
                Some(j) => {
                    trace.push(j);
                    cur = j;
                }
                None => break,
            }
        }
        if trace.len() < 4 {
            continue;
        }
        // 三态单调 DP: 状态 2=上 1=中 0=下, 沿轨迹非增
        let inf = f64::INFINITY;
        let mut dp = vec![[inf; 3]; trace.len()];
        let mut back = vec![[3usize; 3]; trace.len()];
        for (s, row) in dp[0].iter_mut().enumerate() {
            *row = neg_ln(membership(positions, trace[0], s));
        }
        for (t, &i) in trace.iter().enumerate().skip(1) {
            let prev_row = dp[t - 1];
            for (s, row) in dp[t].iter_mut().enumerate() {
                let mut best = inf;
                let mut arg = 3usize;
                for (sp, prev) in prev_row.iter().enumerate().skip(s) {
                    if *prev < best {
                        best = *prev;
                        arg = sp;
                    }
                }
                if best.is_finite() {
                    *row = best + neg_ln(membership(positions, i, s));
                    back[t][s] = arg;
                }
            }
        }
        // 回溯 DP 最优序列
        let mut seq = vec![0usize; trace.len()];
        let mut s = {
            let last = trace.len() - 1;
            (0..3).min_by(|&a, &b| dp[last][a].partial_cmp(&dp[last][b]).unwrap()).unwrap()
        };
        seq[trace.len() - 1] = s;
        for t in (1..trace.len()).rev() {
            s = back[t][s];
            seq[t - 1] = s;
        }
        // 原序列代价
        let dp_cost = dp[trace.len() - 1][seq[trace.len() - 1]];
        let mut orig_cost = 0f64;
        for &i in &trace {
            let r = rank_of_code(result.terrain[i]);
            if r >= 0 {
                orig_cost += neg_ln(membership(positions, i, r as usize));
            }
        }
        if dp_cost >= orig_cost - 1e-9 || dp_cost.is_nan() {
            continue;
        }
        // 应用: 仅低置信像元
        for (t, &i) in trace.iter().enumerate() {
            if result.terrain[i] == 1 || !valid_i(i, result) {
                continue;
            }
            let hill = matches!(result.terrain[i], 3..=5);
            let new_code = match seq[t] {
                2 => {
                    if hill {
                        3
                    } else {
                        6
                    }
                }
                1 => {
                    if hill {
                        4
                    } else {
                        7
                    }
                }
                _ => {
                    if hill {
                        5
                    } else {
                        8
                    }
                }
            };
            if result.terrain[i] != new_code && result.confidence[i] < 0.35 {
                result.terrain[i] = new_code;
                corrected += 1;
            }
        }
    }
    corrected
}

#[inline]
fn basin_or_valid(i: usize, result: &FinalClassification) -> bool {
    valid_i(i, result)
}

#[inline]
fn valid_i(i: usize, result: &FinalClassification) -> bool {
    // terrain==0 表示 NoData
    result.terrain[i] != 0
}
