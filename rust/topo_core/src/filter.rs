//! 并行分离滤波: 均值 / 最小 / 最大（nearest 外推语义, 与 scipy 一致）

use rayon::prelude::*;

/// 跨线程裸指针包装: 各并行迭代项写互不相交的行, 安全性由调用布局保证
#[derive(Clone, Copy)]
struct SendPtr<P>(P);
unsafe impl<P> Send for SendPtr<P> {}
unsafe impl<P> Sync for SendPtr<P> {}
impl<P> std::ops::Deref for SendPtr<P> {
    type Target = P;
    fn deref(&self) -> &P {
        &self.0
    }
}

// 行内写入辅助: along_rows 时 index = l*w+i, 否则 = i*w+l（行间不相交, 并行安全）
#[inline]
unsafe fn put(dst: *mut f32, l: usize, i: usize, w: usize, along_rows: bool, v: f32) {
    let idx = if along_rows { l * w + i } else { i * w + l };
    *dst.add(idx) = v;
}

#[inline]
unsafe fn get(src: *const f32, l: usize, i: isize, w: usize, line_len: usize, along_rows: bool) -> f32 {
    let ii = i.clamp(0, line_len as isize - 1) as usize;
    let idx = if along_rows { l * w + ii } else { ii * w + l };
    *src.add(idx)
}

/// 一维滑动均值（nearest 外推: 窗口恒为 k 宽, 越界以边缘值重复计）
fn mean_axis(src: &[f32], dst: &mut [f32], h: usize, w: usize, k: usize, along_rows: bool) {
    let r = (k / 2) as isize;
    let line_len = if along_rows { w } else { h };
    let n_lines = if along_rows { h } else { w };
    let src_p = SendPtr(src.as_ptr());
    let dst_p = SendPtr(dst.as_mut_ptr());
    (0..n_lines).into_par_iter().for_each(|l| {
        let src_p = *src_p;
        let dst_p = *dst_p;
        let mut ps = vec![0f64; line_len + 1];
        for i in 0..line_len {
            ps[i + 1] =
                ps[i] + unsafe { get(src_p, l, i as isize, w, line_len, along_rows) } as f64;
        }
        let edge_l = ps[1];
        let edge_r = ps[line_len] - ps[line_len - 1];
        for i in 0..line_len as isize {
            let a = i - r;
            let b = i + r;
            let acl = a.clamp(0, line_len as isize - 1) as usize;
            let bcl = b.clamp(0, line_len as isize - 1) as usize;
            let mut sum = ps[bcl + 1] - ps[acl];
            if a < 0 {
                sum += edge_l * (-a) as f64;
            }
            if b >= line_len as isize {
                sum += edge_r * (b - line_len as isize + 1) as f64;
            }
            unsafe {
                put(dst_p, l, i as usize, w, along_rows, (sum / (b - a + 1) as f64) as f32);
            }
        }
    });
}

/// 一维滑动 min/max（nearest 外推: 越界以边缘值参与比较, min/max 与重复计数无关）
fn minmax_axis(
    src: &[f32],
    dst: &mut [f32],
    h: usize,
    w: usize,
    k: usize,
    along_rows: bool,
    want_max: bool,
) {
    let r = (k / 2) as isize;
    let line_len = if along_rows { w } else { h };
    let n_lines = if along_rows { h } else { w };
    let src_p = SendPtr(src.as_ptr());
    let dst_p = SendPtr(dst.as_mut_ptr());
    let cmp = |a: f32, b: f32| if want_max { a > b } else { a < b };
    (0..n_lines).into_par_iter().for_each(|l| {
        let src_p = *src_p;
        let dst_p = *dst_p;
        let g = |i: isize| -> f32 { unsafe { get(src_p, l, i, w, line_len, along_rows) } };
        for i in 0..line_len as isize {
            // 窗口 [i-r, i+r], get 对越界索引做边缘外推
            let a = i - r;
            let b = i + r;
            let acl = a.clamp(0, line_len as isize - 1);
            let bcl = b.clamp(0, line_len as isize - 1);
            let mut m = g(acl);
            for j in (acl + 1)..=bcl {
                let v = g(j);
                if cmp(v, m) {
                    m = v;
                }
            }
            // 窗口越界部分以边缘值重复: min/max 已包含在 clamp 后区间端点中
            if a < 0 {
                let v = g(a);
                if cmp(v, m) {
                    m = v;
                }
            }
            if b >= line_len as isize {
                let v = g(b);
                if cmp(v, m) {
                    m = v;
                }
            }
            unsafe { put(dst_p, l, i as usize, w, along_rows, m) };
        }
    });
}

/// 矩形窗口均值滤波（分离两轴）
pub fn focal_mean(src: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    let k = win | 1;
    let mut tmp = vec![0f32; w * h];
    let mut dst = vec![0f32; w * h];
    mean_axis(src, &mut tmp, h, w, k, true);
    mean_axis(&tmp, &mut dst, h, w, k, false);
    dst
}

/// 矩形窗口最小值滤波
pub fn focal_min(src: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    let k = win | 1;
    let mut tmp = vec![0f32; w * h];
    let mut dst = vec![0f32; w * h];
    minmax_axis(src, &mut tmp, h, w, k, true, false);
    minmax_axis(&tmp, &mut dst, h, w, k, false, false);
    dst
}

/// 矩形窗口最大值滤波
pub fn focal_max(src: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    let k = win | 1;
    let mut tmp = vec![0f32; w * h];
    let mut dst = vec![0f32; w * h];
    minmax_axis(src, &mut tmp, h, w, k, true, true);
    minmax_axis(&tmp, &mut dst, h, w, k, false, true);
    dst
}

/// 局部起伏度 = max - min
pub fn focal_relief(src: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    let mn = focal_min(src, w, h, win);
    let mx = focal_max(src, w, h, win);
    mx.iter().zip(mn.iter()).map(|(a, b)| a - b).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_nearest_center() {
        let w = 4usize;
        let h = 3usize;
        let src: Vec<f32> = (0..w * h).map(|i| i as f32).collect();
        let out = focal_mean(&src, w, h, 3);
        // 中心 (1,1): 3x3 = {0,1,2,4,5,6,8,9,10} 均值 = 5.0
        assert!((out[1 * w + 1] - 5.0).abs() < 1e-4);
        // 角 (0,0): nearest 外推 3x3 = {0,0,1,0,0,1,4,4,5} 均值 = 15/9
        assert!((out[0] - 15.0 / 9.0).abs() < 1e-4);
    }

    #[test]
    fn minmax_nearest() {
        let w = 5usize;
        let h = 5usize;
        let mut src = vec![0f32; w * h];
        src[2 * w + 2] = 9.0;
        let mx = focal_max(&src, w, h, 3);
        assert_eq!(mx[2 * w + 2], 9.0);
        assert_eq!(mx[0], 0.0);
        let mn = focal_min(&src, w, h, 3);
        assert_eq!(mn[2 * w + 2], 0.0);
    }

    #[test]
    fn relief() {
        let w = 6usize;
        let h = 6usize;
        let mut src = vec![1f32; w * h];
        src[3 * w + 3] = 50.0;
        src[0] = -20.0;
        let r = focal_relief(&src, w, h, 3);
        assert_eq!(r[3 * w + 3], 49.0);
    }
}

/// 焦点标准差(方窗, 对齐 ArcGIS FocalStatistics STD; f64 累计防精度损失。
/// O(n·win²), 仅用于粗层小窗(如 TPI 坡位的 101@5m 等效 21@25m))
pub fn focal_std(src: &[f32], w: usize, h: usize, win: usize) -> Vec<f32> {
    let r = (win / 2) as isize;
    let n_sq = (win * win) as f64;
    let mut out = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let mut sum = 0f64;
            let mut sum_sq = 0f64;
            for dy in -r..=r {
                let yy = (y as isize + dy).clamp(0, h as isize - 1) as usize;
                for dx in -r..=r {
                    let xx = (x as isize + dx).clamp(0, w as isize - 1) as usize;
                    let v = src[yy * w + xx] as f64;
                    sum += v;
                    sum_sq += v * v;
                }
            }
            let mean = sum / n_sq;
            let var = (sum_sq / n_sq - mean * mean).max(0.0);
            out[y * w + x] = var.sqrt() as f32;
        }
    }
    out
}
