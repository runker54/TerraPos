//! 精确欧氏距离变换 (EDT) + 最近源索引
//! Meijster (2004) 两遍扫描: 列向扫描得各像元到同列最近源的垂直距离 G,
//! 行向整数 Sep 包络得精确平方距离; 最近源索引由包络中心列回查该列最近源得到。
//! 纯整数运算, 与 scipy.ndimage.distance_transform_edt(..., return_indices=True) 一致。

/// 无源列的占位距离上限, 保证平方运算不溢出 (真实列内距离 ≤ 图高)
const GCAP: i64 = 1 << 20;

#[inline]
fn floordiv(a: i64, b: i64) -> i64 {
    let q = a / b;
    if (a % b != 0) && ((a < 0) != (b < 0)) {
        q - 1
    } else {
        q
    }
}

/// mask=true 处为源。返回 (每像元最近源的线性索引, 欧氏距离)。
pub fn edt_with_index(mask: &[bool], w: usize, h: usize) -> (Vec<u32>, Vec<f32>) {
    let n = w * h;
    let mut best_idx = vec![0u32; n];
    let mut sq = vec![GCAP * GCAP; n];

    // 全图无源: 距离无穷大
    if !mask.iter().any(|b| *b) {
        let dist = vec![f32::INFINITY; n];
        return (best_idx, dist);
    }

    // ---------- Pass 1: 列向, g[i]=到同列最近源的垂直距离, col_src[i]=该源线性索引 ----------
    let mut g = vec![GCAP; n];
    let mut col_src = vec![0u32; n];
    for x in 0..w {
        let mut d: i64 = GCAP;
        let mut src: i64 = -1;
        for y in 0..h {
            let i = y * w + x;
            if mask[i] {
                d = 0;
                src = i as i64;
            }
            g[i] = d;
            if src >= 0 {
                col_src[i] = src as u32;
            }
            if d < GCAP {
                d += 1;
            }
        }
        d = GCAP;
        src = -1;
        for y in (0..h).rev() {
            let i = y * w + x;
            if mask[i] {
                d = 0;
                src = i as i64;
            }
            if d < g[i] {
                g[i] = d;
                if src >= 0 {
                    col_src[i] = src as u32;
                }
            }
            if d < GCAP {
                d += 1;
            }
        }
    }

    // ---------- Pass 2: 行向 Meijster Sep 包络 ----------
    let mut v = vec![0usize; w];
    let mut z = vec![0i64; w + 1];
    for y in 0..h {
        let row = y * w;
        let f: Vec<i64> = (0..w).map(|x| g[row + x] * g[row + x]).collect();
        let mut q: usize = 0;
        v[0] = 0;
        z[0] = i64::MIN;
        z[1] = i64::MAX;
        for x in 1..w {
            let xi = x as i64;
            let vi = v[q] as i64;
            let mut s = floordiv(
                (xi * xi + f[x]) - (vi * vi + f[v[q]]),
                2 * xi - 2 * vi,
            );
            while s <= z[q] {
                q -= 1;
                let vi2 = v[q] as i64;
                s = floordiv(
                    (xi * xi + f[x]) - (vi2 * vi2 + f[v[q]]),
                    2 * xi - 2 * vi2,
                );
            }
            q += 1;
            v[q] = x;
            z[q] = s;
            z[q + 1] = i64::MAX;
        }
        let mut qi: usize = 0;
        for x in 0..w {
            while z[qi + 1] < x as i64 {
                qi += 1;
            }
            let c = v[qi];
            let dx = x as i64 - c as i64;
            sq[row + x] = dx * dx + f[c];
            best_idx[row + x] = col_src[row + c];
        }
    }

    let dist: Vec<f32> = sq.iter().map(|s| (*s as f64).sqrt() as f32).collect();
    (best_idx, dist)
}

/// 用最近有效像元值填充 nodata (非有限值)
pub fn fill_nodata_nearest(dem: &mut [f32], w: usize, h: usize) {
    let invalid: Vec<bool> = dem.iter().map(|v| !v.is_finite()).collect();
    if !invalid.iter().any(|b| *b) {
        return;
    }
    let valid: Vec<bool> = invalid.iter().map(|b| !b).collect();
    let (vidx, _) = edt_with_index(&valid, w, h);
    for (i, bad) in invalid.iter().enumerate() {
        if *bad {
            dem[i] = dem[vidx[i] as usize];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edt_single_source() {
        let w = 5usize;
        let h = 5usize;
        let mut mask = vec![false; w * h];
        mask[2 * w + 2] = true;
        let (idx, dist) = edt_with_index(&mask, w, h);
        assert_eq!(idx[2 * w + 2], (2 * w + 2) as u32);
        assert_eq!(dist[2 * w + 2], 0.0);
        assert_eq!(dist[2 * w + 4], 2.0);
        assert!((dist[0] - 8f32.sqrt()).abs() < 1e-5);
        assert_eq!(idx[0], (2 * w + 2) as u32);
    }

    #[test]
    fn edt_two_sources() {
        let w = 5usize;
        let h = 1usize;
        let mut mask = vec![false; w];
        mask[0] = true;
        mask[4] = true;
        let (idx, dist) = edt_with_index(&mask, w, h);
        assert_eq!(dist[2], 2.0);
        assert_eq!(dist[1], 1.0);
        assert_eq!(dist[3], 1.0);
        assert_eq!(idx[1], 0);
        assert_eq!(idx[3], 4);
    }

    #[test]
    fn edt_matches_bruteforce() {
        // 与暴力最近邻逐像元比对 (5x7 含分散源)
        let w = 5usize;
        let h = 7usize;
        let mut mask = vec![false; w * h];
        let sources = [3usize, 9, 17, 30];
        for &s in &sources {
            mask[s] = true;
        }
        let (idx, dist) = edt_with_index(&mask, w, h);
        for y in 0..h {
            for x in 0..w {
                let i = y * w + x;
                let dsums: Vec<(f64, u32)> = sources
                    .iter()
                    .map(|&s| {
                        let (sx, sy) = (s % w, s / w);
                        (
                            ((x as f64 - sx as f64).powi(2) + (y as f64 - sy as f64).powi(2))
                                .sqrt(),
                            s as u32,
                        )
                    })
                    .collect();
                let bd = dsums.iter().map(|(d, _)| *d).fold(f64::INFINITY, f64::min);
                // 距离必须精确一致
                assert!((dist[i] - bd as f32).abs() < 1e-4, "dist mismatch at {i}");
                // 平局时任意一个最近源都合法
                assert!(
                    dsums
                        .iter()
                        .any(|&(d, s)| d - bd < 1e-9 && s == idx[i]),
                    "idx mismatch at ({x},{y}): got {}, allowed {:?}",
                    idx[i],
                    dsums
                        .iter()
                        .filter(|&&(d, _)| d - bd < 1e-9)
                        .map(|&(_, s)| s)
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    #[test]
    fn edt_big_image_sanity() {
        // 1.4 亿像元级大图: 单源中心, 远处像元距离必须精确
        let (w, h) = (20000usize, 7000usize);
        let mut mask = vec![false; w * h];
        let (cx, cy) = (w / 2, h / 2);
        mask[cy * w + cx] = true;
        let (idx, dist) = edt_with_index(&mask, w, h);
        // 探针: 源右侧 100 像元
        let probe = cy * w + cx + 100;
        assert!((dist[probe] - 100.0).abs() < 1.5, "dist[probe]={}", dist[probe]);
        // 探针: 源左下角方向
        let probe2 = (cy - 50) * w + cx - 50;
        let expect = ((50.0f64).powi(2) * 2.0).sqrt() as f32;
        assert!((dist[probe2] - expect).abs() < 1.5, "dist[probe2]={}", dist[probe2]);
        assert_eq!(idx[probe], (cy * w + cx) as u32);
    }

    #[test]
    fn edt_no_source_degenerate() {
        let w = 3usize;
        let h = 3usize;
        let mask = vec![false; w * h];
        let (idx, dist) = edt_with_index(&mask, w, h);
        assert!(dist.iter().all(|d| d.is_infinite()));
        assert!(idx.iter().all(|i| *i == 0));
    }
}
