//! 多尺度地貌因子金字塔: 稳健 DEV / 起伏 H / NDEV 与特征尺度选择。
//!
//! 所有尺度以米定义, 在消费边界(`usable_scales`/滤波内部)一次换算为像元。
//! 只存粗层网格上的逐尺度层, 避免在原生分辨率保存多套大数组。

use crate::error::{CoreError, Result};
use crate::filter::focal_robust_stats_valid;
use crate::input::RasterShape;

/// 候选尺度族(米)
pub const BASE_SCALES_M: [f64; 6] = [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0];

/// 单尺度地貌因子层
#[derive(Debug, Clone)]
pub struct ScaleLayer {
    pub radius_m: f64,
    /// DEV = z - local_median
    pub deviation_m: Vec<f32>,
    /// H = P95 - P05(米)
    pub relief_m: Vec<f32>,
    /// NDEV = DEV / (1.4826 * MAD + 0.1)
    pub normalized_deviation: Vec<f32>,
}

/// 多尺度金字塔与特征尺度选择结果
#[derive(Debug, Clone)]
pub struct ScalePyramid {
    pub layers: Vec<ScaleLayer>,
    /// 每像元特征分析尺度(米)
    pub characteristic_scale_m: Vec<f32>,
    /// 每像元特征尺度对应的稳健起伏(米)
    pub characteristic_relief_m: Vec<f32>,
    pub shape: RasterShape,
}

/// 由网格范围裁剪候选尺度: 删除小于 5 像元与大于有效区短边四分之一的尺度
pub fn usable_scales(shape: RasterShape) -> Vec<f64> {
    let min_scale = 5.0 * shape.resolution_m;
    let short_side = shape.width.min(shape.height) as f64 * shape.resolution_m;
    let max_scale = short_side * 0.25;
    BASE_SCALES_M
        .iter()
        .copied()
        .filter(|&s| s >= min_scale && s <= max_scale)
        .collect()
}

/// 构建尺度金字塔并选择特征尺度。
/// 选择规则: 首个连续两个增长率 g < growth_threshold 且 H > 5 m 的尺度;
/// 无满足者时在 H > 5 m 的层中取 |g - threshold| 最小; 仍无则取最小可用尺度。
pub fn build_scale_pyramid(
    dem: &[f32],
    valid: &[bool],
    shape: RasterShape,
    growth_threshold: f32,
) -> Result<ScalePyramid> {
    let scales = usable_scales(shape);
    if scales.len() < 3 {
        return Err(CoreError::Invalid(format!(
            "DEM extent too small for multiscale classification: 只剩 {} 个可用尺度 {:?}",
            scales.len(),
            scales
        )));
    }
    let n = shape.width * shape.height;
    let mut reliefs: Vec<Vec<f32>> = Vec::with_capacity(scales.len());
    let mut layers: Vec<ScaleLayer> = Vec::with_capacity(scales.len());
    for &radius in &scales {
        let (med, mad, p05, p95) = focal_robust_stats_valid(dem, valid, shape, radius);
        let mut dev = vec![0f32; n];
        let mut ndev = vec![0f32; n];
        let mut relief = vec![0f32; n];
        for i in 0..n {
            if !valid[i] || !med[i].is_finite() {
                continue;
            }
            let d = dem[i] - med[i];
            dev[i] = d;
            ndev[i] = d / (1.4826 * mad[i] + 0.1);
            relief[i] = p95[i] - p05[i];
        }
        reliefs.push(relief.clone());
        layers.push(ScaleLayer {
            radius_m: radius,
            deviation_m: dev,
            relief_m: relief,
            normalized_deviation: ndev,
        });
    }
    // 尺度相邻增长率 g_j = (H_{j+1}-H_j)/max(H_{j+1}, 0.1)
    let ng = scales.len() - 1;
    let mut growth = Vec::with_capacity(ng);
    for j in 0..ng {
        let (h0, h1) = (&reliefs[j], &reliefs[j + 1]);
        let mut g = vec![0f32; n];
        for i in 0..n {
            let (a, b) = (h0[i], h1[i]);
            g[i] = (b - a) / b.max(0.1);
        }
        growth.push(g);
    }

    let mut characteristic_scale_m = vec![0f32; n];
    let mut characteristic_relief_m = vec![0f32; n];
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        // 首个连续两个 transition 收敛且 H>5m 的尺度
        let mut chosen: Option<usize> = None;
        for j in 0..ng.saturating_sub(1) {
            if growth[j][i] < growth_threshold
                && growth[j + 1][i] < growth_threshold
                && reliefs[j][i] > 5.0
            {
                chosen = Some(j);
                break;
            }
        }
        // 回退: H>5m 的层中 |g-th| 最小(g 序列比层数短 1, 末层用最后 g)
        if chosen.is_none() {
            let mut best: Option<(usize, f32)> = None;
            for j in 0..ng {
                if reliefs[j][i] > 5.0 {
                    let score = (growth[j][i] - growth_threshold).abs();
                    if best.map(|(_, bs)| score < bs).unwrap_or(true) {
                        best = Some((j, score));
                    }
                }
            }
            chosen = best.map(|(j, _)| j);
        }
        // 最终回退: 最小可用尺度
        let j = chosen.unwrap_or(0);
        characteristic_scale_m[i] = scales[j] as f32;
        characteristic_relief_m[i] = reliefs[j][i];
    }

    Ok(ScalePyramid {
        layers,
        characteristic_scale_m,
        characteristic_relief_m,
        shape,
    })
}

/// 地貌亚类尺度上下限(米)。口径与项目现有规则一致:
/// <500 m 丘陵(2000 m 窗口起伏 <200 m 为低丘), 500-1000 低山,
/// 1000-3500 中山, 3500-5000 高山, >=5000 极高山。
pub fn scale_bounds_for(elevation_m: f32, relief_m: f32, data_max_m: f32) -> (f32, f32) {
    const HILL_MAX: f32 = 500.0;
    const LOW_HILL_RELIEF: f32 = 200.0;
    if elevation_m < HILL_MAX {
        if relief_m < LOW_HILL_RELIEF {
            (125.0, 750.0)
        } else {
            (250.0, 1500.0)
        }
    } else if elevation_m < 1000.0 {
        (250.0, 2000.0)
    } else if elevation_m < 3500.0 {
        (500.0, 4000.0)
    } else if elevation_m < 5000.0 {
        (750.0, 6000.0)
    } else {
        (1000.0, 8000.0_f32.min(data_max_m).max(1000.0))
    }
}
