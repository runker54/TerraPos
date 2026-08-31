//! 地貌个体分割: 峰顶提取(非极大抑制) + marker watershed + 连通域标记

use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

const N8: [(isize, isize); 8] = [
    (-1, -1), (0, -1), (1, -1),
    (-1, 0), (1, 0),
    (-1, 1), (0, 1), (1, 1),
];

/// 将候选中的等值平台收缩为代表元(平台质心最近像元, 平局取索引小者)。
/// 局部极大平台的所有像元都是候选; 收缩到质心使种子落在峰顶中心。
fn plateau_centers(smooth: &[f32], w: usize, cands: &[(f32, usize)]) -> Vec<(f32, usize)> {
    let mut visited = vec![false; smooth.len()];
    let mut is_cand = vec![false; smooth.len()];
    for &(_, i) in cands {
        is_cand[i] = true;
    }
    let mut out = Vec::new();
    for &(_, i) in cands {
        if visited[i] {
            continue;
        }
        let val = smooth[i];
        let mut comp = vec![i];
        visited[i] = true;
        let mut qi = 0;
        while qi < comp.len() {
            let c = comp[qi];
            qi += 1;
            let cx = (c % w) as isize;
            let cy = (c / w) as isize;
            for (dx, dy) in N8 {
                let nx = cx + dx;
                let ny = cy + dy;
                if nx < 0 || ny < 0 || nx >= w as isize || ny >= (smooth.len() / w) as isize {
                    continue;
                }
                let j = ny as usize * w + nx as usize;
                if !visited[j] && is_cand[j] && smooth[j] == val {
                    visited[j] = true;
                    comp.push(j);
                }
            }
        }
        let n = comp.len() as f64;
        let gx = comp.iter().map(|&c| (c % w) as f64).sum::<f64>() / n;
        let gy = comp.iter().map(|&c| (c / w) as f64).sum::<f64>() / n;
        let best = *comp
            .iter()
            .min_by(|&&a, &&b| {
                let da = ((a % w) as f64 - gx).powi(2) + ((a / w) as f64 - gy).powi(2);
                let db = ((b % w) as f64 - gx).powi(2) + ((b / w) as f64 - gy).powi(2);
                da.partial_cmp(&db).unwrap().then(a.cmp(&b))
            })
            .unwrap();
        out.push((val, best));
    }
    out
}

/// 格网地形突起度 (topographic prominence): 峰顶到连通更高地的最高鞍点的高差。
/// 降序并查集: 格点按高程降序加入, 与已处理 8 邻域的集合合并; 集合首次被他集
/// 合并处的当前点高程即其最高鞍点。突出度天然多尺度——大山突出度大、小山突出度
/// 小, 无需窗口参数; 返回 (每像元突出度, 是否峰顶/平台代表)。
/// 全图最高点的突出度以相对全图最低点计。
pub fn prominence_map(dem: &[f32], w: usize, h: usize) -> (Vec<f32>, Vec<bool>) {
    let n = w * h;
    let zmin = dem.iter().copied().fold(f32::INFINITY, f32::min);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| dem[b].partial_cmp(&dem[a]).unwrap().then(a.cmp(&b)));

    let mut parent: Vec<u32> = (0..n as u32).collect();
    let mut processed = vec![false; n];
    let mut peak_z = vec![0f32; n]; // 集合节点: 集合内最高点高程
    let mut prom = vec![f32::NAN; n]; // 集合节点: 已定突出度(NAN=未定)
    let mut is_peak = vec![false; n];
    let mut node_of = vec![u32::MAX; n]; // 峰顶 -> 其集合创建节点(突出度所在)

    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        loop {
            let p = parent[x as usize];
            if p == x {
                return x;
            }
            let g = parent[p as usize];
            parent[x as usize] = g;
            x = p;
        }
    }

    for &p in &order {
        let zp = dem[p];
        // 已处理邻居的集合(去重)
        let px = (p % w) as isize;
        let py = (p / w) as isize;
        let mut adj = [0u32; 8];
        let mut na = 0usize;
        for (dx, dy) in N8 {
            let nx = px + dx;
            let ny = py + dy;
            if nx < 0 || ny < 0 || nx >= w as isize || ny >= h as isize {
                continue;
            }
            let q = ny as usize * w + nx as usize;
            if processed[q] {
                let r = find(&mut parent, q as u32);
                if !adj[..na].contains(&r) {
                    adj[na] = r;
                    na += 1;
                }
            }
        }
        if na == 0 {
            is_peak[p] = true;
            peak_z[p] = zp;
            node_of[p] = p as u32;
        } else {
            // 主集合 = 峰顶最高者, 其余集合在此割点并入(记录突出度)
            let mut rmax = adj[0];
            for &r in &adj[1..na] {
                if peak_z[r as usize] > peak_z[rmax as usize] {
                    rmax = r;
                }
            }
            parent[p] = rmax;
            for &r in &adj[..na] {
                if r != rmax {
                    if prom[r as usize].is_nan() {
                        prom[r as usize] = peak_z[r as usize] - zp;
                    }
                    parent[r as usize] = rmax;
                }
            }
        }
        processed[p] = true;
    }

    // 突出度定在峰顶创建节点上; 未被他集合合并(全图最高)者取相对全图最低点
    let mut prom_map = vec![0f32; n];
    for i in 0..n {
        if is_peak[i] {
            let r = node_of[i] as usize;
            let v = prom[r];
            prom_map[i] = if v.is_nan() { peak_z[r] - zmin } else { v };
        }
    }
    (prom_map, is_peak)
}

/// 突出度种子: 突出度 >= `min_prom_m` 的峰顶/平台代表, 无需窗口与 NMS。
pub fn prominence_seeds(smooth: &[f32], w: usize, h: usize, min_prom_m: f32) -> Vec<usize> {
    let (prom, is_peak) = prominence_map(smooth, w, h);
    (0..w * h)
        .into_par_iter()
        .filter(|&i| is_peak[i] && prom[i] >= min_prom_m)
        .map(|i| i)
        .collect()
}

/// 多尺度 TPI 特征尺度投票(尺度空间语义, 无种子无窗口选择):
/// 每像元在尺度族 `scales_cells`(像元)上计算 TPI, 取 TPI 最大的尺度为其特征尺度;
/// 最大 TPI > 0 的像元视为山体, 8 连通域直接作为地貌个体单元。
/// 返回每像元单元标签(0 = 非山体/平原)。
pub fn multiscale_units(dem: &[f32], w: usize, h: usize, scales_cells: &[usize]) -> Vec<u32> {
    let n = w * h;
    let mut best_tpi = vec![f32::NEG_INFINITY; n];
    let mut best_scale = vec![0u8; n];
    for (si, &win) in scales_cells.iter().enumerate() {
        let win = win.max(1) | 1;
        let t = crate::terrain::tpi(dem, w, h, win);
        for i in 0..n {
            if t[i] > best_tpi[i] {
                best_tpi[i] = t[i];
                best_scale[i] = si as u8;
            }
        }
    }
    let mask: Vec<bool> = best_tpi.iter().map(|&t| t > 0.0).collect();
    let (lab, _) = label_connected(&mask, w, h);
    lab.iter().map(|&l| l.max(0) as u32).collect()
}

/// 混合种子: 突出度语义 ∪ 距离语义。
/// 突出度 ≥ min_prom_m 的峰直接保留(干净个体);
/// NMS 窗口种子中, 距所有突出度种子 ≥ 半窗者补入(保留坡面鼓包等
/// 突出度不足但距离上独立的个体)。两语义互补, 兼顾干净与完整。
pub fn hybrid_seeds(
    smooth: &[f32],
    w: usize,
    h: usize,
    min_dist_cells: usize,
    min_prom_m: f32,
) -> Vec<usize> {
    let prom_seeds = prominence_seeds(smooth, w, h, min_prom_m);
    let fixed = peak_seeds(smooth, w, h, min_dist_cells);
    let half2 = (min_dist_cells / 2).max(1);
    let half2 = (half2 * half2) as f64;
    let mut seeds = prom_seeds.clone();
    for &f in &fixed {
        let (fx, fy) = ((f % w) as f64, (f / w) as f64);
        let near_prom = prom_seeds.iter().any(|&p| {
            let (px, py) = ((p % w) as f64, (p / w) as f64);
            (fx - px) * (fx - px) + (fy - py) * (fy - py) < half2
        });
        if !near_prom {
            seeds.push(f);
        }
    }
    seeds
}



/// 峰顶种子提取: 平滑面 DEM 的局部极大 + 非极大抑制(最小间距=窗口)
/// 语义与 skimage.feature.peak_local_max(min_distance=...) 一致:
/// 候选 = focal_max(2*min_dist+1) 的相等点; 按高程降序贪心, 抑制欧氏距离 < min_dist 的邻点。
pub fn peak_seeds(
    smooth: &[f32],
    w: usize,
    h: usize,
    min_dist_cells: usize,
) -> Vec<usize> {
    let win = (min_dist_cells * 2 + 1) | 1;
    let mx = crate::filter::focal_max(smooth, w, h, win);
    let mn = crate::filter::focal_min(smooth, w, h, win);
    // 候选: 局部极大且窗口内严格突出(排除零突出度的平坦基准面)
    let mut cands: Vec<(f32, usize)> = (0..w * h)
        .into_par_iter()
        .filter(|&i| smooth[i] >= mx[i] && smooth[i] > mn[i])
        .map(|i| (smooth[i], i))
        .collect();
    // 高程降序、索引升序（平顶取先者）
    cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
    let mut cands = plateau_centers(smooth, w, &cands);
    cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
    let md2 = (min_dist_cells * min_dist_cells) as f64;
    let mut chosen: Vec<(f64, f64)> = Vec::new();
    let mut seeds = Vec::new();
    for (_, i) in cands {
        let x = (i % w) as f64;
        let y = (i / w) as f64;
        let near = chosen
            .iter()
            .any(|(cx, cy)| (x - cx) * (x - cx) + (y - cy) * (y - cy) < md2);
        if !near {
            chosen.push((x, y));
            seeds.push(i);
        }
    }
    seeds
}

/// 分亚类多尺度峰顶提取: 每个亚类用其"地貌个体尺度"窗口做局部极大+NMS,
/// 不同亚类区域独立提取(窗口随地貌空间自适应), 合并为全域峰顶种子。
/// `windows`: (亚类编码, 峰顶最小间距像元) 查表。
pub fn peak_seeds_zoned(
    smooth: &[f32],
    w: usize,
    h: usize,
    sub: &[u8],
    windows: &[(u8, usize)],
) -> Vec<usize> {
    let mut all = Vec::new();
    for (sub_code, r_cells) in windows {
        // 该亚类区域置顶, 区域外 -inf(不参与局部极大)
        let mut zone = smooth.to_vec();
        for (i, &s) in sub.iter().enumerate() {
            if s != *sub_code {
                zone[i] = f32::NEG_INFINITY;
            }
        }
        let win = (r_cells * 2 + 1) | 1;
        // 区域外 -inf 供 focal_max; +inf 供 focal_min(突出度只在区域内度量)
        let mut zone_hi = vec![f32::INFINITY; w * h];
        for (i, &s) in sub.iter().enumerate() {
            if s == *sub_code {
                zone_hi[i] = smooth[i];
            }
        }
        let mx = crate::filter::focal_max(&zone, w, h, win);
        let mn = crate::filter::focal_min(&zone_hi, w, h, win);
        // 候选: 区域内局部极大且窗口内严格突出(排除零突出度的平坦基准面), 排除 -inf 占位
        let mut cands: Vec<(f32, usize)> = (0..w * h)
            .into_par_iter()
            .filter(|&i| zone[i] != f32::NEG_INFINITY && zone[i] >= mx[i] && zone[i] > mn[i])
            .map(|i| (zone[i], i))
            .collect();
        // 高程降序贪心 NMS(欧氏距离 >= 窗口); 等值平台先收缩到质心代表元
        cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
        let mut cands = plateau_centers(smooth, w, &cands);
        cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(&b.1)));
        let md2 = (*r_cells * *r_cells) as f64;
        let mut chosen: Vec<(f64, f64)> = Vec::new();
        for (_, i) in cands {
            let x = (i % w) as f64;
            let y = (i / w) as f64;
            if chosen
                .iter()
                .all(|(cx, cy)| (x - cx) * (x - cx) + (y - cy) * (y - cy) >= md2)
            {
                chosen.push((x, y));
                all.push(i);
            }
        }
    }
    all
}

/// marker watershed: 从峰顶种子沿负地形从高到低淹没(等价 skimage.segmentation.watershed(-dem, markers))
/// 返回每像元的单元标签(0 = 无)
pub fn watershed(dem: &[f32], w: usize, h: usize, seeds: &[usize]) -> Vec<u32> {
    let n = w * h;
    let mut labels = vec![0u32; n];
    let mut heap: BinaryHeap<Reverse<(i64, u32, u32)>> = BinaryHeap::with_capacity(1024);
    let ordered = crate::hydro::ordered;
    for (si, &s) in seeds.iter().enumerate() {
        labels[s] = (si + 1) as u32;
        heap.push(Reverse((ordered(-dem[s]), (si + 1) as u32, s as u32)));
    }
    while let Some(Reverse((_e, lab, i))) = heap.pop() {
        let xi = (i as usize % w) as isize;
        let yi = (i as usize / w) as isize;
        for (dx, dy) in N8 {
            let nx = xi + dx;
            let ny = yi + dy;
            if nx < 0 || ny < 0 || nx >= w as isize || ny >= h as isize {
                continue;
            }
            let j = (ny as usize) * w + nx as usize;
            if labels[j] == 0 {
                labels[j] = lab;
                heap.push(Reverse((ordered(-dem[j]), lab, j as u32)));
            }
        }
    }
    labels
}

/// 连通域标记（8 邻域, BFS, 标签 1..n; 0=背景）
pub fn label_connected(mask: &[bool], w: usize, h: usize) -> (Vec<i32>, usize) {
    let n = w * h;
    let mut lab = vec![0i32; n];
    let mut cur = 0i32;
    let mut stack: Vec<usize> = Vec::with_capacity(1024);
    for s in 0..n {
        if !mask[s] || lab[s] != 0 {
            continue;
        }
        cur += 1;
        lab[s] = cur;
        stack.push(s);
        while let Some(i) = stack.pop() {
            let x = (i % w) as isize;
            let y = (i / w) as isize;
            for (dx, dy) in N8 {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= w as isize || ny >= h as isize {
                    continue;
                }
                let j = (ny as usize) * w + nx as usize;
                if mask[j] && lab[j] == 0 {
                    lab[j] = cur;
                    stack.push(j);
                }
            }
        }
    }
    (lab, cur as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn watershed_markers_all_labeled() {
        // 6x1: 种子 idx2(峰20) 与 idx5(谷10); 淹没后全像元有归属
        let w = 6usize;
        let h = 1usize;
        let dem = vec![10.0f32, 5.0, 20.0, 10.0, 5.0, 10.0];
        let seeds = vec![2usize, 5usize];
        let labels = watershed(&dem, w, h, &seeds);
        // 种子自身标签
        assert_eq!(labels[2], 1);
        assert_eq!(labels[5], 2);
        // 全域被标记
        assert!(labels.iter().all(|&l| l > 0));
        // 峰顶 20 归种子 idx2; 右端谷底(idx5 邻域) 归种子 idx5
        assert_eq!(labels[3], 1); // 峰 20 的直接邻居归峰种子
        assert_eq!(labels[5], 2);
    }

    #[test]
    fn prominence_two_peaks_hierarchy() {
        // 9x1: 主峰(4)=100, 鞍(5)=80, 次峰(6)=90, 其余 50
        // 次峰通往更高地的最低门槛 = 鞍 80 -> 突出度 = 90-80 = 10
        // 主峰为全图最高点 -> 突出度 = 100 - 最低点(50) = 50
        let w = 9usize;
        let h = 1usize;
        let mut dem = vec![50.0f32; w];
        dem[4] = 100.0;
        dem[5] = 80.0;
        dem[6] = 90.0;
        let (prom, is_peak) = prominence_map(&dem, w, h);
        assert!(is_peak[4] && is_peak[6]);
        assert!(!is_peak[3] && !is_peak[5]);
        assert!((prom[6] - 10.0).abs() < 1e-4, "prom[6]={}", prom[6]);
        assert!((prom[4] - 50.0).abs() < 1e-4, "prom[4]={}", prom[4]);
        // 阈值 15 只留主峰; 阈值 5 两峰都留
        assert_eq!(prominence_seeds(&dem, w, h, 15.0), vec![4]);
        let mut s = prominence_seeds(&dem, w, h, 5.0);
        s.sort_unstable();
        assert_eq!(s, vec![4, 6]);
    }

    #[test]
    fn prominence_equal_peaks_share_saddle() {
        // 等高双峰共享鞍 50: 双方突出度 = 100-50 = 50
        let w = 3usize;
        let h = 1usize;
        let dem = vec![100.0f32, 50.0, 100.0];
        let (prom, is_peak) = prominence_map(&dem, w, h);
        assert!(is_peak[0] && is_peak[2]);
        assert!((prom[0] - 50.0).abs() < 1e-4, "prom[0]={}", prom[0]);
        assert!((prom[2] - 50.0).abs() < 1e-4, "prom[2]={}", prom[2]);
    }

    #[test]
    fn prominence_flat_base_no_seeds() {
        // 全平: 唯一平台代表突出度 = peak - zmin = 0 -> 被阈值滤除
        let w = 4usize;
        let h = 4usize;
        let dem = vec![7.0f32; w * h];
        let seeds = prominence_seeds(&dem, w, h, 5.0);
        assert!(seeds.is_empty());
    }

    #[test]
    fn prominence_hill_on_slope() {
        // 斜坡 100-x 上的小丘(9..11 抬升 20): 丘顶 idx9=111
        // 丘的突出度 = 111 - 丘脚门槛(上坡侧 92); 坡下端点(0)突出度 = 100-92 = 8
        let w = 20usize;
        let h = 1usize;
        let mut dem: Vec<f32> = (0..w).map(|x| 100.0 - x as f32).collect();
        for x in [9usize, 10, 11] {
            dem[x] += 20.0;
        }
        let (prom, is_peak) = prominence_map(&dem, w, h);
        assert!(is_peak[9]);
        // 丘顶 111 为全图最高点 -> 突出度相对全图最低点(81) = 30
        assert!((prom[9] - 30.0).abs() < 1e-4, "prom[9]={}", prom[9]);
        // 坡足端点(0)=100 经丘脚门槛 92 与更高地合并 -> 突出度 8
        assert!((prom[0] - 8.0).abs() < 1e-4, "prom[0]={}", prom[0]);
        // 阈值 10: 只有丘顶; 边界坡足(prom=8)被滤除
        assert_eq!(prominence_seeds(&dem, w, h, 10.0), vec![9]);
    }

    #[test]
    fn multiscale_units_two_hills() {
        // 两个不同尺度的山丘 + 平原: 正 TPI 连通域各自成单元, 平原为 0
        let w = 60usize;
        let h = 30usize;
        let mut dem = vec![100.0f32; w * h];
        let gauss = |x: f64, y: f64, cx: f64, cy: f64, a: f64, s: f64| {
            a * (-((x - cx) * (x - cx) + (y - cy) * (y - cy)) / (2.0 * s * s)).exp()
        };
        for y in 0..h {
            for x in 0..w {
                let (fx, fy) = (x as f64, y as f64);
                let z = gauss(fx, fy, 15.0, 15.0, 40.0, 6.0)
                    + gauss(fx, fy, 45.0, 15.0, 40.0, 14.0);
                dem[y * w + x] += z as f32;
            }
        }
        let units = multiscale_units(&dem, w, h, &[5, 11, 21]);
        assert_eq!(units[0], 0, "平原应为背景");
        assert!(units[15 * w + 15] > 0, "小丘中心应为山体");
        assert!(units[15 * w + 45] > 0, "大丘中心应为山体");
        // 两丘之间的鞍部/远处平原不应与丘同单元
        assert_ne!(units[15 * w + 15], 0);
    }

    #[test]
    fn label_8connected() {
        let w = 3usize;
        let h = 3usize;
        let mut mask = vec![false; w * h];
        mask[0] = true;
        mask[4] = true; // (1,1) 对角 -> 8 邻域与 (0,0) 连通
        let (lab, n) = label_connected(&mask, w, h);
        assert_eq!(n, 1);
        assert_eq!(lab[0], 1);
        assert_eq!(lab[4], 1);
    }

    #[test]
    fn peak_seeds_two_peaks() {
        // 8x1: 峰 (1)=10, (5)=9, 谷 3
        let w = 8usize;
        let h = 1usize;
        let mut dem = vec![3.0f32; w];
        dem[1] = 10.0;
        dem[5] = 9.0;
        let seeds = peak_seeds(&dem, w, h, 2);
        assert!(seeds.contains(&1));
        assert!(seeds.contains(&5));
        assert_eq!(seeds.len(), 2);
    }

    #[test]
    fn peak_seeds_zoned_three_basins() {
        // h=1 全低山: 三峰 (10,30) (55,26) (100,28), 峰间 45 像元 > 方窗全宽 41
        // 各峰独立成种子(验证 windows 查表与 zone 提取主流程)
        use crate::filter::focal_mean;
        let w = 140usize;
        let h = 1usize;
        let mut dem = vec![2.0f32; w];
        dem[10] = 30.0;
        dem[55] = 26.0;
        dem[100] = 28.0;
        let sub = vec![3u8; w];
        let smooth = focal_mean(&dem, w, h, 3);
        let windows = [(3u8, 20usize)];
        let seeds = peak_seeds_zoned(&smooth, w, h, &sub, &windows);
        // 三峰各自成种子(索引即列号)
        assert!(seeds.contains(&10), "seeds={:?}", seeds);
        assert!(seeds.contains(&55), "seeds={:?}", seeds);
        assert!(seeds.contains(&100), "seeds={:?}", seeds);
    }
}
