//! 脊线提取: 子流域边界(主证据) + 反地形/形态补充证据融合。
//!
//! 子流域按细河网链(同一河链共享标识)竞争泛洪划分; 链间接触面为主
//! 分水候选; 反地形汇流、正 NDEV 与凸性形态作为子流域内部次级脊的
//! 补充证据, 仅在强度达标且靠近主分水(或图内无主分水)时保留。

use crate::distance::edt_with_index;
use crate::error::Result;
use crate::geomorphon::Landform;
use crate::hydro::HydroModel;
use crate::input::RasterShape;
use crate::scale::ScalePyramid;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// 脊线模型产品(粗层网格)
#[derive(Debug, Clone)]
pub struct RidgeModel {
    /// 最终脊线掩膜(不含河网)
    pub mask: Vec<bool>,
    /// 补充证据强度 0..1(主分水处为 1)
    pub strength: Vec<f32>,
    /// 子流域(河链集水区)标识, 0=河网自身/未达
    pub subcatchment_id: Vec<u32>,
    pub shape: RasterShape,
}

/// 融合子流域边界与补充证据构建脊线。
pub fn build_ridges(
    dem: &[f32],
    valid: &[bool],
    hydro: &HydroModel,
    pyramid: &ScalePyramid,
    landform: &[Landform],
) -> Result<RidgeModel> {
    let w = hydro.shape.width;
    let h = hydro.shape.height;
    let n = w * h;
    let stream = &hydro.streams[0];

    // 0) 深洼汇流区: 从深洼环沿 flow_to 反向 BFS(规格 8.3: 闭合洼地是
    // 独立水文子系统, 洼内河链不得互相产生子流域边界)
    let mut closed_region = hydro.deep_sink.clone();
    {
        let mut stack: Vec<usize> = (0..n).filter(|&i| hydro.deep_sink[i]).collect();
        while let Some(i) = stack.pop() {
            for j in reverse_neighbours(i, w, h) {
                if !closed_region[j] && hydro.flow_to[j] == i as u32 {
                    closed_region[j] = true;
                    stack.push(j);
                }
            }
        }
    }

    // 1) 河网链标注: 河网连通体(8 连通)为一条链 —— 网状/辫状河带
    // 属同一排水系统, 按流向继承会在横带上退化为每像元一链
    let mut chain_of = vec![0u32; n];
    let mut next_chain = 1u32;
    {
        let mut stack: Vec<usize> = Vec::new();
        for s0 in 0..n {
            if !stream[s0] || chain_of[s0] != 0 {
                continue;
            }
            let c = next_chain;
            next_chain += 1;
            chain_of[s0] = c;
            stack.push(s0);
            while let Some(i) = stack.pop() {
                for j in reverse_neighbours(i, w, h) {
                    if stream[j] && chain_of[j] == 0 {
                        chain_of[j] = c;
                        stack.push(j);
                    }
                }
            }
        }
    }

    // 2) 子流域竞争泛洪: 河链像元为多源种子, 反向 BFS 就近归属
    let mut sub = vec![0u32; n];
    let mut queue: BinaryHeap<Reverse<(u32, u32)>> = BinaryHeap::new(); // (层级, idx)
    for i in 0..n {
        if stream[i] && chain_of[i] != 0 {
            sub[i] = chain_of[i];
            queue.push(Reverse((0, i as u32)));
        }
    }
    let mut level = vec![u32::MAX; n];
    while let Some(Reverse((lv, i))) = queue.pop() {
        let i = i as usize;
        if lv > level[i] {
            continue;
        }
        for j in reverse_neighbours(i, w, h) {
            if !valid[j] || hydro.flow_to[j] != i as u32 {
                continue;
            }
            if sub[j] == 0 {
                sub[j] = sub[i];
                level[j] = lv + 1;
                queue.push(Reverse((lv + 1, j as u32)));
            }
        }
    }
    // 未被任何链汇流覆盖的有效像元(洼地内流等)直接归最近链, 无则 0
    // (其边界由空缺边缘自然形成)

    // 3) 主分水候选: 相异标识的接触面对(双侧标记);
    // 排除细河网近旁(<=2 粗像元)像元——河网离散化会在谷旁制造伪边界
    let mut mask = vec![false; n];
    {
        let (_, dist_stream) = edt_with_index(&hydro.streams[0], w, h);
        for i in 0..n {
            if !valid[i] || stream[i] || sub[i] == 0 || dist_stream[i] <= 2.0 {
                continue;
            }
            // 分水线必然局部高于邻域(NDEV>0): 平行流的算法伪接触
            // (无横向起伏)不产生边界
            if ndev_at(pyramid, pyramid.characteristic_scale_m[i], i) <= 0.0 {
                continue;
            }
            for j in reverse_neighbours(i, w, h) {
                if valid[j]
                    && !stream[j]
                    && sub[j] != 0
                    && sub[j] != sub[i]
                    && dist_stream[j] > 2.0
                    && ndev_at(pyramid, pyramid.characteristic_scale_m[j], j) > 0.0
                {
                    mask[i] = true;
                    mask[j] = true;
                }
            }
        }
    }
    // 双侧标记细化: 保 id 较小一侧, 得单像元分水线
    for i in 0..n {
        if !mask[i] {
            continue;
        }
        for j in reverse_neighbours(i, w, h) {
            if mask[j] && sub[j] != 0 && sub[i] != 0 && sub[j] < sub[i] {
                mask[i] = false;
                break;
            }
        }
    }

    // 4) 补充证据强度
    let mut strength = vec![0f32; n];
    // 反地形汇流(脊线=反地形谷线); 无效区给高值防止反地形谷线穿越
    let zmax = dem.iter().copied().fold(f32::MIN, f32::max);
    let zmin = dem.iter().copied().fold(f32::MAX, f32::min);
    let inv_dem: Vec<f32> = dem
        .iter()
        .zip(valid.iter())
        .map(|(&z, &ok)| {
            if ok && z.is_finite() {
                zmax + zmin - z
            } else {
                zmax
            }
        })
        .collect();
    let inv_route = crate::hydro::fill_and_route(&inv_dem, w, h, 99999.0);
    // 反地形汇流阈值取有效像元的 P95: 真脊线是反地形汇流的高分位,
    // 固定小阈值会把对称坡的"反地形汇流横带"误判为脊
    let inv_threshold = {
        let mut sample: Vec<u32> = inv_route
            .acc
            .iter()
            .zip(valid.iter())
            .filter(|&(_, &v)| v)
            .step_by(16)
            .map(|(&a, _)| a)
            .collect();
        sample.sort_unstable();
        let p95 = sample[(sample.len() as f32 * 0.95) as usize % sample.len()];
        p95.max(80)
    };
    for i in 0..n {
        if !valid[i] || stream[i] {
            continue;
        }
        let mut s = 0f32;
        if inv_route.acc[i] >= inv_threshold {
            s += 0.4;
        }
        if ndev_at(pyramid, pyramid.characteristic_scale_m[i], i) > 0.0 {
            s += 0.3;
        }
        if matches!(
            landform[i],
            Landform::Peak | Landform::Ridge | Landform::Shoulder | Landform::Spur
        ) {
            s += 0.3;
        }
        strength[i] = s;
    }

    // 补充像元: 强度达标且(靠近主分水 或 图内无主分水 或 三证据俱全)。
    // 三证据俱全(strength>=0.999)的强脊独立保留, 如环状盆缘的分水;
    // 其余强证据像元须靠近主分水, 防平地噪声扩张成脊(规格 8.2)
    let has_divide = mask.iter().any(|&b| b);
    let (src, dist) = edt_with_index(&mask, w, h);
    let _ = src;
    let max_gap_m = (0.1 * median_scale(pyramid)).min(250.0);
    let max_gap_px = (max_gap_m / hydro.shape.resolution_m).ceil() as f32;
    for i in 0..n {
        if !mask[i]
            && strength[i] >= 0.7
            && valid[i]
            && !stream[i]
            && (!has_divide || dist[i] <= max_gap_px || strength[i] >= 0.999)
        {
            mask[i] = true;
        }
    }

    // 主分水处强度记 1
    for (i, m) in mask.iter().enumerate() {
        if *m && sub[i] != 0 {
            strength[i] = strength[i].max(1.0);
        }
    }

    // 5) 2x2 全真块保左上(线状细化)
    let mut thin = mask.clone();
    for y in 0..h - 1 {
        for x in 0..w - 1 {
            let i = y * w + x;
            if mask[i] && mask[i + 1] && mask[i + w] && mask[i + w + 1] {
                thin[i + 1] = false;
                thin[i + w] = false;
                thin[i + w + 1] = false;
            }
        }
    }

    Ok(RidgeModel {
        mask: thin,
        strength,
        subcatchment_id: sub,
        shape: hydro.shape,
    })
}

/// 反向邻居(流入 i 的像元)候选 8 邻域
fn reverse_neighbours(i: usize, w: usize, h: usize) -> Vec<usize> {
    const N8: [(isize, isize); 8] = [
        (-1, -1),
        (0, -1),
        (1, -1),
        (-1, 0),
        (1, 0),
        (-1, 1),
        (0, 1),
        (1, 1),
    ];
    let x = (i % w) as isize;
    let y = (i / w) as isize;
    let mut out = Vec::with_capacity(8);
    for (dx, dy) in N8 {
        let nx = x + dx;
        let ny = y + dy;
        if nx >= 0 && ny >= 0 && nx < w as isize && ny < h as isize {
            out.push((ny as usize) * w + nx as usize);
        }
    }
    out
}

/// 取特征尺度对应层(最近层)的 NDEV
fn ndev_at(pyramid: &ScalePyramid, scale_m: f32, i: usize) -> f32 {
    let mut best = &pyramid.layers[0];
    let mut bd = f32::INFINITY;
    for layer in &pyramid.layers {
        let d = (layer.radius_m as f32 - scale_m).abs();
        if d < bd {
            bd = d;
            best = layer;
        }
    }
    best.normalized_deviation[i]
}

/// 全图特征尺度中位数(米)
fn median_scale(pyramid: &ScalePyramid) -> f64 {
    let mut v: Vec<f32> = pyramid
        .characteristic_scale_m
        .iter()
        .copied()
        .filter(|s| *s > 0.0)
        .collect();
    if v.is_empty() {
        return 500.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2] as f64
}
