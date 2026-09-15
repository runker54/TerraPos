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
        shape: hydro.shape,
    })
}
