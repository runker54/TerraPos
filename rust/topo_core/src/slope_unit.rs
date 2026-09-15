//! 坡面单元标注: 河网与脊线双屏障 + 坡向相容泛洪。
//!
//! 单元由"细河网(谷屏障) + 脊线(脊屏障)"围合, 同单元内坡向圆差
//! 不超过 90 度; 过小单元并入共享边界最长的相邻单元。

use crate::error::Result;
use crate::hydro::HydroModel;
use crate::input::RasterShape;
use crate::ridge::RidgeModel;
use std::collections::HashMap;

/// 坡面单元产品(粗层网格)
#[derive(Debug, Clone)]
pub struct SlopeUnits {
    /// 单元标识(1 起, 0=屏障/无效/未分配)
    pub unit_id: Vec<u32>,
    /// 下坡向 8 扇区(0..7), 8=平坡向不敏感
    pub aspect_sector: Vec<u8>,
    /// 谷屏障(细河网)
    pub valley_mask: Vec<bool>,
    /// 脊屏障
    pub ridge_mask: Vec<bool>,
    /// 河链集水区标识(透传自脊线模型, 供谷底对象聚合)
    pub subcatchment_id: Vec<u32>,
    pub shape: RasterShape,
}

/// 构建坡面单元。
pub fn build_slope_units(
    dem: &[f32],
    valid: &[bool],
    hydro: &HydroModel,
    ridges: &RidgeModel,
) -> Result<SlopeUnits> {
    let w = hydro.shape.width;
    let h = hydro.shape.height;
    let n = w * h;
    let valley_mask = hydro.streams[0].clone();
    let ridge_mask = ridges.mask.clone();

    // 坡向扇区: 中心差分下坡方向, 45 度扇区; 平坡向不敏感(8)
    let res = hydro.shape.resolution_m;
    let mut aspect_sector = vec![8u8; n];
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let i = y * w + x;
            if !valid[i] || valley_mask[i] || ridge_mask[i] {
                continue;
            }
            let dzdx = (dem[i + 1] - dem[i - 1]) as f64 / (2.0 * res);
            let dzdy = (dem[i + w] - dem[i - w]) as f64 / (2.0 * res);
            let mag = dzdx.hypot(dzdy);
            if mag < 1e-6 {
                continue; // 平地: 扇区 8
            }
            // 下坡方向 = -梯度方向; 扇区 0=东, 逆时针每 45 度
            let ang = (-dzdy).atan2(-dzdx).to_degrees();
            let norm = (ang + 360.0).rem_euclid(360.0);
            let sector = ((norm / 45.0).round() as i64).rem_euclid(8) as u8;
            aspect_sector[i] = sector;
        }
    }

    // 相容泛洪: 坡向圆差 <= 90 度(<=2 扇区), 平坡向(8)任意相容
    let compat = |a: u8, b: u8| -> bool {
        if a == 8 || b == 8 {
            return true;
        }
        let d = (a as i64 - b as i64).abs().min(8 - (a as i64 - b as i64).abs());
        d <= 2
    };
    let mut unit_id = vec![0u32; n];
    let mut cur = 0u32;
    let mut sizes: HashMap<u32, usize> = HashMap::new();
    for s0 in 0..n {
        if !valid[s0] || valley_mask[s0] || ridge_mask[s0] || unit_id[s0] != 0 {
            continue;
        }
        cur += 1;
        unit_id[s0] = cur;
        let mut stack = vec![s0];
        let mut sz = 0usize;
        while let Some(i) = stack.pop() {
            sz += 1;
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
                if !valid[j]
                    || valley_mask[j]
                    || ridge_mask[j]
                    || unit_id[j] != 0
                {
                    continue;
                }
                // 子流域约束: 不同子流域(分水两侧)不得进入同一单元,
                // 脊线掩膜的局部缺口不再导致跨脊单元
                if ridges.subcatchment_id[i] != 0
                    && ridges.subcatchment_id[j] != 0
                    && ridges.subcatchment_id[i] != ridges.subcatchment_id[j]
                {
                    continue;
                }
                if compat(aspect_sector[i], aspect_sector[j]) {
                    unit_id[j] = cur;
                    stack.push(j);
                }
            }
        }
        sizes.insert(cur, sz);
    }

    // 过小单元(<4 粗像元)并入共享边界最长的相邻单元
    let mut merged = sizes.clone();
    let small: Vec<u32> = sizes
        .iter()
        .filter(|(&_, &s)| s < 4)
        .map(|(&id, _)| id)
        .collect();
    for id in small {
        // 统计与各相邻单元的共享边界数
        let mut contact: HashMap<u32, usize> = HashMap::new();
        for i in 0..n {
            if unit_id[i] != id {
                continue;
            }
            let x = (i % w) as i64;
            let y = (i / w) as i64;
            for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1)] {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                let other = unit_id[j];
                if other != 0 && other != id && merged.get(&other).copied().unwrap_or(0) >= 4 {
                    *contact.entry(other).or_default() += 1;
                }
            }
        }
        if let Some((&target, _)) = contact.iter().max_by_key(|(&t, &c)| (c, t)) {
            for v in unit_id.iter_mut() {
                if *v == id {
                    *v = target;
                }
            }
            let sz = merged.remove(&id).unwrap_or(0);
            *merged.entry(target).or_insert(0) += sz;
        }
        // 无合格邻居时保留小单元(总比无归属好)
    }

    Ok(SlopeUnits {
        unit_id,
        aspect_sector,
        valley_mask,
        ridge_mask,
        subcatchment_id: ridges.subcatchment_id.clone(),
        shape: hydro.shape,
    })
}

// ---------------- 约束距离/自适应尺度/相对位置(Task 7) ----------------

use crate::hydro;
use crate::scale::{scale_bounds_for, ScalePyramid};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// 坡面几何产品(粗层网格)
#[derive(Debug, Clone)]
pub struct SlopeGeometry {
    /// 谷/脊锚点像元索引(u32::MAX=缺失)
    pub valley_anchor: Vec<u32>,
    pub ridge_anchor: Vec<u32>,
    /// 约束距离(米, 不跨单元屏障)
    pub distance_to_valley_m: Vec<f32>,
    pub distance_to_ridge_m: Vec<f32>,
    /// 局部坡面宽 L = dv + dr
    pub local_width_m: Vec<f32>,
    /// 最终分析尺度 R*(米)
    pub adaptive_scale_m: Vec<f32>,
    /// 流向关联 HAND(米)
    pub hand_m: Vec<f32>,
    /// qd / qz / 融合位置量 q
    pub q_distance: Vec<f32>,
    pub q_elevation: Vec<f32>,
    pub relative_position: Vec<f32>,
    /// 锚点缺失或单元不闭合(消费方据此降低置信度)
    pub low_confidence: Vec<bool>,
    pub shape: RasterShape,
}

/// 屏障感知多源 Dijkstra: 种子向外传播, 不同单元的邻居不入队;
/// 边代价 = Δs * (1 + 2*|sin(无向坡向差)|), 沿坡向移动代价最低。
fn constrained_dijkstra(
    seeds: &[bool],
    units: &SlopeUnits,
    valid: &[bool],
) -> (Vec<f32>, Vec<u32>) {
    let w = units.shape.width;
    let h = units.shape.height;
    let res = units.shape.resolution_m;
    let n = w * h;
    let mut dist = vec![f32::INFINITY; n];
    let mut anchor = vec![u32::MAX; n];
    let mut heap: BinaryHeap<Reverse<(i64, u32)>> = BinaryHeap::with_capacity(1024);
    for i in 0..n {
        if seeds[i] {
            dist[i] = 0.0;
            anchor[i] = i as u32;
            heap.push(Reverse((hydro::ordered(0.0), i as u32)));
        }
    }
    // 移动方向扇区(与 aspect_sector 同约定: 0=东, 每 45 度)
    let move_sector = |dx: i64, dy: i64| -> u8 {
        let ang = (-dy as f64).atan2(dx as f64).to_degrees();
        (((ang + 360.0).rem_euclid(360.0) / 45.0).round() as i64).rem_euclid(8) as u8
    };
    while let Some(Reverse((dbits, i))) = heap.pop() {
        let i = i as usize;
        if hydro::ordered(dist[i]) > dbits {
            continue; // lazy deletion
        }
        // 屏障像元(脊/谷, unit 0)只作锚定终点, 不得作为中转,
        // 否则约束距离会跨脊/跨谷把对侧单元连起来
        let i_is_barrier = units.unit_id[i] == 0;
        if i_is_barrier && !seeds[i] {
            continue;
        }
        let x = (i % w) as i64;
        let y = (i / w) as i64;
        for (dx, dy) in [(-1i64, 0i64), (1, 0), (0, -1), (0, 1), (-1, -1), (1, -1), (-1, 1), (1, 1)] {
            let nx = x + dx;
            let ny = y + dy;
            if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                continue;
            }
            let j = (ny as usize) * w + nx as usize;
            if !valid[j] {
                continue;
            }
            // 跨单元禁止: 两侧均为单元像元且不同(屏障由"不入队"阻断)
            if !i_is_barrier
                && units.unit_id[j] != 0
                && units.unit_id[i] != units.unit_id[j]
            {
                continue;
            }
            let step = res * if dx != 0 && dy != 0 { std::f64::consts::SQRT_2 } else { 1.0 };
            let diff = {
                let a = units.aspect_sector[i];
                let b = move_sector(dx, dy);
                if a == 8 {
                    0.0
                } else {
                    let d = (a as i64 - b as i64).abs().min(8 - (a as i64 - b as i64).abs());
                    (d as f64) * 45.0f64.to_radians()
                }
            };
            let nd = dist[i] + (step * (1.0 + 2.0 * diff.sin().abs())) as f32;
            if nd < dist[j] {
                dist[j] = nd;
                anchor[j] = anchor[i];
                if units.unit_id[j] != 0 {
                    heap.push(Reverse((hydro::ordered(nd), j as u32)));
                }
            }
        }
    }
    (dist, anchor)
}

/// 计算约束几何: 逐河网等级 Dijkstra -> L_k -> 等级选择 -> R* -> 相对位置 q。
pub fn build_slope_geometry(
    dem: &[f32],
    valid: &[bool],
    units: &SlopeUnits,
    hydro: &hydro::HydroModel,
    pyramid: &ScalePyramid,
    low_relief_m: f32,
) -> Result<SlopeGeometry> {
    let n = units.shape.width * units.shape.height;
    // 脊距离与锚(一次)
    let (dr, ranchor) = constrained_dijkstra(&units.ridge_mask, units, valid);
    // 逐等级谷距离与锚
    let mut dv_k: [Vec<f32>; 4] = Default::default();
    let mut va_k: [Vec<u32>; 4] = Default::default();
    for k in 0..4 {
        let (d, a) = constrained_dijkstra(&hydro.streams[k], units, valid);
        dv_k[k] = d;
        va_k[k] = a;
    }

    let data_max = dem.iter().copied().fold(f32::MIN, f32::max);
    let mut geom = SlopeGeometry {
        valley_anchor: va_k[0].clone(),
        ridge_anchor: ranchor,
        distance_to_valley_m: dv_k[0].clone(),
        distance_to_ridge_m: dr.clone(),
        local_width_m: vec![f32::NAN; n],
        adaptive_scale_m: vec![f32::NAN; n],
        hand_m: vec![f32::NAN; n],
        q_distance: vec![f32::NAN; n],
        q_elevation: vec![f32::NAN; n],
        relative_position: vec![f32::NAN; n],
        low_confidence: vec![false; n],
        shape: units.shape,
    };
    // HAND(粗层)
    if let Ok(hand) = hydro::hand_to_stream(&hydro.conditioned, &hydro.flow_to, &hydro.streams[0], &vec![true; n]) {
        geom.hand_m = hand;
    }

    for i in 0..n {
        if !valid[i] {
            continue;
        }
        let dr_i = dr[i];
        // 等级选择: min |ln(L_k/(2R0))| + 0.5|ln(L_{k+1}/L_k)| (末级末项为 0)
        let r0 = pyramid.characteristic_scale_m[i].max(1.0) as f64;
        let mut best_k = 0usize;
        let mut best_score = f64::INFINITY;
        for k in 0..4 {
            let lv = dv_k[k][i];
            if !lv.is_finite() {
                continue;
            }
            let l_k = (lv + dr_i) as f64;
            if l_k <= 0.0 {
                continue;
            }
            let mut score = (l_k / (2.0 * r0)).ln().abs();
            if k + 1 < 4 && dv_k[k + 1][i].is_finite() {
                let l_next = (dv_k[k + 1][i] + dr_i) as f64;
                if l_next > 0.0 {
                    score += 0.5 * (l_next / l_k).ln().abs();
                }
            }
            if score < best_score {
                best_score = score;
                best_k = k;
            }
        }
        // 所选等级的 dv/锚(全部缺失时回退细等级 k=0), 诊断层与 L/q 保持一致
        let (dv_i, va_i) = if dv_k[best_k][i].is_finite() {
            (dv_k[best_k][i], va_k[best_k][i])
        } else {
            (dv_k[0][i], va_k[0][i])
        };
        geom.distance_to_valley_m[i] = dv_i;
        geom.valley_anchor[i] = va_i;
        let l = if dv_i.is_finite() && dr_i.is_finite() {
            dv_i + dr_i
        } else {
            f32::NAN
        };
        geom.local_width_m[i] = l;
        // R* = clip(0.5 L, Rmin, Rmax)
        let (rmin, rmax) = scale_bounds_for(
            dem[i],
            pyramid.characteristic_relief_m[i],
            data_max,
        );
        geom.adaptive_scale_m[i] = if l.is_finite() {
            (0.5 * l).clamp(rmin, rmax)
        } else {
            rmax
        };
        // 锚点缺失 -> 低置信
        let has_v = va_i != u32::MAX;
        let has_r = geom.ridge_anchor[i] != u32::MAX;
        geom.low_confidence[i] = !(has_v && has_r);
        // qd / qz / q
        let zv = if has_v { dem[va_i as usize] } else { f32::NAN };
        let zr = if has_r {
            dem[geom.ridge_anchor[i] as usize]
        } else {
            f32::NAN
        };
        let qd = if has_v && l.is_finite() {
            (dv_i / (dv_i + dr_i + 0.001)).clamp(0.0, 1.0)
        } else {
            f32::NAN
        };
        let qz = if has_v && has_r {
            ((dem[i] - zv) / (zr - zv + 0.001)).clamp(0.0, 1.0)
        } else {
            f32::NAN
        };
        geom.q_distance[i] = qd;
        geom.q_elevation[i] = qz;
        let h_star = pyramid.characteristic_relief_m[i];
        let w_distance = 0.45f32
            + (0.75 - 0.45) * ((low_relief_m - h_star) / low_relief_m).clamp(0.0, 1.0);
        let q = match (qd.is_finite(), qz.is_finite()) {
            (true, true) => (w_distance * qd + (1.0 - w_distance) * qz).clamp(0.0, 1.0),
            (true, false) => qd,
            (false, true) => qz,
            (false, false) => f32::NAN,
        };
        geom.relative_position[i] = q;
    }
    Ok(geom)
}
