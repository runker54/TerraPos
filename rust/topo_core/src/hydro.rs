//! 水文分析: Priority-Flood 填洼(带 z-limit) + D8 流向 + 汇流累积
//!
//! 算法: Barnes et al. 2014 "Priority-Flood" 变体。
//! - 填洼与流向一步完成: 每个像元首次从堆中弹出的邻居即其下游(flat 安全);
//! - z-limit 语义与 ESRI Fill 一致: 洼地深度<=z_limit 完全填平, 超限完全不填
//!   (深洼成为内流区);
//! - 汇流累积按弹出顺序逆序累加(天然拓扑序, 无需独立排序)。
//!
//! 性能: BinaryHeap<Reverse<(i64 有序高度, u32 idx)>>, O(N log N);
//! 本区 25m 层 537 万像元实测秒级。

use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(Debug, Clone)]
pub struct FlowResult {
    /// 填洼后表面
    pub filled: Vec<f32>,
    /// 每像元下游像元的一维索引(边界出流指向自身)
    pub flow_to: Vec<u32>,
    /// Priority-Flood 弹出顺序(汇流累积的拓扑序)
    pub pop_order: Vec<u32>,
    /// 超出 z-limit 的深洼像元(内流区)
    pub is_deep_sink: Vec<bool>,
    /// 汇流累积像元数(含自身)
    pub acc: Vec<u32>,
}

const N8: [(isize, isize); 8] = [
    (-1, -1), (0, -1), (1, -1),
    (-1, 0), (1, 0),
    (-1, 1), (0, 1), (1, 1),
];

/// f32 -> 可全序比较的 i64(IEEE 位变换; 数据经 nodata 填充保证有限)
#[inline]
pub fn ordered(v: f32) -> i64 {
    let b = v.to_bits() as i32;
    if b >= 0 { b as i64 } else { !(b as i64) }
}

/// 填洼 + D8 流向 + 汇流累积
pub fn fill_and_route(dem: &[f32], w: usize, h: usize, z_limit: f32) -> FlowResult {
    let n = w * h;
    let mut filled = dem.to_vec();
    let mut flow_to = vec![u32::MAX; n];
    let mut visited = vec![false; n];
    let mut is_deep_sink = vec![false; n];
    let mut pop_order: Vec<u32> = Vec::with_capacity(n);
    // 最小堆: Reverse((有序高度 i64, idx u32))
    let mut heap: BinaryHeap<Reverse<(i64, u32)>> = BinaryHeap::with_capacity(1024);

    // 1) 边界像元为初始种子
    for y in 0..h {
        for x in 0..w {
            if x == 0 || y == 0 || x == w - 1 || y == h - 1 {
                let i = y * w + x;
                visited[i] = true;
                flow_to[i] = i as u32;
                heap.push(Reverse((ordered(dem[i]), i as u32)));
            }
        }
    }

    // 2) Priority-Flood: 低处先弹出, 邻居首次被发现时确定其下游
    while let Some(Reverse((_e_ord, i))) = heap.pop() {
        let i = i as usize;
        let e = filled[i]; // 该像元的填充面高程
        pop_order.push(i as u32);
        let xi = (i % w) as isize;
        let yi = (i / w) as isize;
        for (dx, dy) in N8 {
            let nx = xi + dx;
            let ny = yi + dy;
            if nx < 0 || ny < 0 || nx >= w as isize || ny >= h as isize {
                continue;
            }
            let j = (ny as usize) * w + nx as usize;
            if visited[j] {
                continue;
            }
            visited[j] = true;
            let fill_e = if e > dem[j] { e } else { dem[j] };
            flow_to[j] = i as u32;
            if fill_e - dem[j] > z_limit {
                // 深洼边缘: 不填充, 以原值入堆(内流区局部海平面)
                is_deep_sink[j] = true;
                filled[j] = dem[j];
                heap.push(Reverse((ordered(dem[j]), j as u32)));
            } else {
                filled[j] = fill_e;
                heap.push(Reverse((ordered(fill_e), j as u32)));
            }
        }
    }

    // 兜底: 未访问像元(不应存在)流向自身
    for (i, f) in flow_to.iter_mut().enumerate() {
        if *f == u32::MAX {
            *f = i as u32;
        }
    }

    // 3) 汇流累积: 逆拓扑序累加(子先于父)
    let mut acc = vec![1u32; n];
    for &i in pop_order.iter().rev() {
        let i = i as usize;
        let t = flow_to[i] as usize;
        if t != i {
            acc[t] += acc[i];
        }
    }

    FlowResult { filled, flow_to, pop_order, is_deep_sink, acc }
}

/// 从汇流累积栅格提取河网掩膜(阈值=集水像元数)
pub fn stream_mask(acc: &[u32], threshold: u32) -> Vec<bool> {
    acc.iter().map(|&a| a >= threshold).collect()
}


/// 验证辅助: 对 flow_to 森林计算每节点的子树大小(含自身)
#[cfg(test)]
fn subtree_size(flow_to: &[u32], n: usize) -> Vec<u32> {
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, &t) in flow_to.iter().enumerate() {
        let t = t as usize;
        if t != i {
            children[t].push(i);
        }
    }
    let mut size = vec![1u32; n];
    // 迭代后序(避免深递归): 用栈
    let mut stack: Vec<(usize, bool)> = Vec::new();
    for root in 0..n {
        if flow_to[root] == root as u32 {
            stack.push((root, false));
            while let Some((node, processed)) = stack.pop() {
                if processed {
                    let mut s = 1u32;
                    for &c in &children[node] {
                        s += size[c];
                    }
                    size[node] = s;
                } else {
                    stack.push((node, true));
                    for &c in &children[node] {
                        stack.push((c, false));
                    }
                }
            }
        }
    }
    size
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fill_shallow_depression() {
        // 4x4 全 10, 中心 2x2 洼地 6(深 4), z_limit=5 -> 完全填平
        let w = 4usize;
        let h = 4usize;
        let mut dem = vec![10f32; w * h];
        for y in 1..3 {
            for x in 1..3 {
                dem[y * w + x] = 6.0;
            }
        }
        let r = fill_and_route(&dem, w, h, 5.0);
        assert!(r.filled.iter().all(|v| (v - 10.0).abs() < 1e-6));
        // acc 语义验证: acc[i] = 以 i 为根的 flow_to 子树大小(含自身)
        assert_eq!(subtree_size(&r.flow_to, w * h), r.acc);
    }

    #[test]
    fn deep_sink_beyond_zlimit() {
        // 4x4 中心 2x2 洼地 2(深 8), z_limit=3 -> 完全不填
        let w = 4usize;
        let h = 4usize;
        let mut dem = vec![10f32; w * h];
        for y in 1..3 {
            for x in 1..3 {
                dem[y * w + x] = 2.0;
            }
        }
        let r = fill_and_route(&dem, w, h, 3.0);
        assert_eq!(r.filled[1 * w + 1], 2.0);
        assert_eq!(r.filled[2 * w + 2], 2.0);
        assert!(r.is_deep_sink.iter().any(|b| *b));
    }

    #[test]
    fn accumulation_mass_balance() {
        // 伪随机地形: 汇流总量守恒 = 像元总数; 流向满足拓扑序
        let w = 16usize;
        let h = 16usize;
        let dem: Vec<f32> = (0..w * h).map(|i| ((i * 7919) % 101) as f32).collect();
        let r = fill_and_route(&dem, w, h, 50.0);
        assert_eq!(subtree_size(&r.flow_to, w * h), r.acc);
        let mut order_pos = vec![0usize; w * h];
        for (p, &i) in r.pop_order.iter().enumerate() {
            order_pos[i as usize] = p;
        }
        for (i, &t) in r.flow_to.iter().enumerate() {
            let t = t as usize;
            if t != i {
                assert!(order_pos[t] < order_pos[i], "流向违反拓扑序 @{}", i);
            }
        }
    }

    #[test]
    fn stream_mask_threshold() {
        let acc = vec![1u32, 5, 100, 3];
        let m = stream_mask(&acc, 5);
        assert_eq!(m, vec![false, true, true, false]);
    }
}
