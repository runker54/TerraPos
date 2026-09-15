//! 输入检查、有效区管理与地貌表面构建。
//!
//! 职责: 校验投影米制元数据、区分外部 NoData 与内部小孔、
//! 填补不超过面积上限的内部孔、生成一次保形平滑的地貌表面 `geomorph`。

use crate::distance::edt_with_index;
use crate::error::{CoreError, Result};
use crate::geotiff::{read_f32, GeoMeta};
use std::path::Path;

/// 栅格形状(分辨率恒为米)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RasterShape {
    pub width: usize,
    pub height: usize,
    pub resolution_m: f64,
}

/// 输入准备配置
#[derive(Debug, Clone)]
pub struct InputConfig {
    /// 允许填补的内部 NoData 孔最大面积(平方米, 0.01 km²)
    pub max_interior_hole_area_m2: f64,
    /// 地貌表面保形平滑半径(米, 不小于 1 像元)
    pub geomorph_smooth_radius_m: f64,
}

impl Default for InputConfig {
    fn default() -> Self {
        InputConfig {
            max_interior_hole_area_m2: 10_000.0,
            geomorph_smooth_radius_m: 5.0,
        }
    }
}

/// 预处理后的 DEM: 原始值(仅接受的内部孔被填)、地貌平滑表面与有效掩膜
#[derive(Debug)]
pub struct PreparedDem {
    pub raw: Vec<f32>,
    pub geomorph: Vec<f32>,
    pub valid: Vec<bool>,
    pub shape: RasterShape,
    pub meta: GeoMeta,
}

/// 校验 GeoTIFF 元数据: 尺寸、正方形正分辨率、投影米制坐标系
pub fn validate_meta(meta: &GeoMeta) -> Result<RasterShape> {
    if meta.width < 3 || meta.height < 3 {
        return Err(CoreError::Invalid(format!(
            "DEM 尺寸过小 ({}x{}), 至少需要 3x3 像元",
            meta.width, meta.height
        )));
    }
    let sx = meta.pixel_scale[0];
    let sy = meta.pixel_scale[1];
    if !(sx > 0.0 && sy > 0.0) {
        return Err(CoreError::Invalid(format!(
            "DEM 像元分辨率必须为正 (观测 {sx} x {sy})"
        )));
    }
    if (sx - sy).abs() / sx.max(sy) > 0.01 {
        return Err(CoreError::Invalid(format!(
            "DEM X/Y 分辨率差异超过 1% (观测 {sx} x {sy})"
        )));
    }
    if !meta.is_projected_metre() {
        let mt = meta
            .geo_key_u16(1024)
            .map(|v| v.to_string())
            .unwrap_or_else(|| "缺失".into());
        let lu = meta
            .geo_key_u16(3076)
            .map(|v| v.to_string())
            .unwrap_or_else(|| "缺失".into());
        return Err(CoreError::Invalid(format!(
            "DEM 必须为投影坐标系且线性单位为米, 请先重投影 \
             (观测 GTModelType={mt}, ProjLinearUnits={lu}, 像元 {sx} m)"
        )));
    }
    Ok(RasterShape {
        width: meta.width as usize,
        height: meta.height as usize,
        resolution_m: sx,
    })
}

/// 由原始高程值构建 PreparedDem:
/// 1. 有限且非 NoData 哨兵的像元进入有效掩膜, 有效比例不足 1% 报错;
/// 2. 与边界连通的无效区永久无效;
/// 3. 内部无效连通对象面积不超过上限时以最近有效值填补(仅此处理改写 raw);
/// 4. 对有效像元做一次 valid-aware 圆形均值平滑得到地貌表面。
pub fn prepare_values(values: Vec<f32>, meta: GeoMeta, cfg: &InputConfig) -> Result<PreparedDem> {
    let shape = validate_meta(&meta)?;
    let n = shape.width * shape.height;
    if values.len() != n {
        return Err(CoreError::Invalid(format!(
            "像元数不匹配: 元数据 {} vs 数据 {}",
            n,
            values.len()
        )));
    }
    let nd = meta.nodata.unwrap_or(f32::NAN);
    let mut valid: Vec<bool> = values.iter().map(|&v| v.is_finite() && v != nd).collect();
    let valid_count = valid.iter().filter(|&&v| v).count();
    if (valid_count as f64) < 0.01 * n as f64 {
        return Err(CoreError::Invalid(format!(
            "DEM 有效像元比例过低: {:.2}% (要求 >= 1%)",
            100.0 * valid_count as f64 / n as f64
        )));
    }
    let mut raw = values;

    // 外部(与图幅边界连通)无效区 flood fill, 其余无效像元记为内部孔
    let interior = mark_interior_holes(&valid, shape);

    // 内部孔按连通对象面积决定是否填补
    let mut fill_idx: Vec<usize> = Vec::new();
    {
        let w = shape.width;
        let mut lab = vec![0u32; n];
        let mut cur = 0u32;
        let mut stack: Vec<usize> = Vec::new();
        for s0 in 0..n {
            if !interior[s0] || lab[s0] != 0 {
                continue;
            }
            cur += 1;
            lab[s0] = cur;
            stack.push(s0);
            let mut comp = vec![s0];
            while let Some(i) = stack.pop() {
                let x = i % w;
                let y = i / w;
                for (dx, dy) in NEIGH4 {
                    let nx = x as i64 + dx;
                    let ny = y as i64 + dy;
                    if nx < 0 || ny < 0 || nx >= w as i64 || ny >= shape.height as i64 {
                        continue;
                    }
                    let j = ny as usize * w + nx as usize;
                    if interior[j] && lab[j] == 0 {
                        lab[j] = cur;
                        stack.push(j);
                        comp.push(j);
                    }
                }
            }
            let area_m2 = comp.len() as f64 * shape.resolution_m * shape.resolution_m;
            if area_m2 <= cfg.max_interior_hole_area_m2 {
                fill_idx.extend(comp);
            }
        }
    }
    if !fill_idx.is_empty() {
        let (src, _) = edt_with_index(&valid, shape.width, shape.height);
        for &i in &fill_idx {
            let s = src[i] as usize;
            raw[i] = raw[s];
            valid[i] = true;
        }
    }

    let geomorph = geomorph_smooth(&raw, &valid, shape, cfg.geomorph_smooth_radius_m);
    Ok(PreparedDem {
        raw,
        geomorph,
        valid,
        shape,
        meta,
    })
}

/// 读取 GeoTIFF 并完成输入准备
pub fn prepare_input(path: &Path, cfg: &InputConfig) -> Result<PreparedDem> {
    let (values, meta) = read_f32(path)?;
    prepare_values(values, meta, cfg)
}

const NEIGH4: [(i64, i64); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];

/// 标记内部孔像元: 无效且不与图幅边界连通
fn mark_interior_holes(valid: &[bool], shape: RasterShape) -> Vec<bool> {
    let w = shape.width;
    let h = shape.height;
    let n = w * h;
    let mut external = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    for x in 0..w {
        for y in [0usize, h - 1] {
            let i = y * w + x;
            if !valid[i] && !external[i] {
                external[i] = true;
                stack.push(i);
            }
        }
    }
    for y in 0..h {
        for x in [0usize, w - 1] {
            let i = y * w + x;
            if !valid[i] && !external[i] {
                external[i] = true;
                stack.push(i);
            }
        }
    }
    while let Some(i) = stack.pop() {
        let x = i % w;
        let y = i / w;
        for (dx, dy) in NEIGH4 {
            let nx = x as i64 + dx;
            let ny = y as i64 + dy;
            if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                continue;
            }
            let j = ny as usize * w + nx as usize;
            if !valid[j] && !external[j] {
                external[j] = true;
                stack.push(j);
            }
        }
    }
    (0..n).map(|i| !valid[i] && !external[i]).collect()
}

/// valid-aware 圆形均值平滑(米制半径 -> 像元半径一次换算)。
/// 无效像元不输出(返回 0, 下游一律以 valid 掩膜过滤); 平面区域保持不变。
fn geomorph_smooth(raw: &[f32], valid: &[bool], shape: RasterShape, radius_m: f64) -> Vec<f32> {
    let r_px = ((radius_m / shape.resolution_m).round() as i64).max(1);
    let r2 = (r_px * r_px) as f64;
    let mut offs: Vec<(i64, i64)> = Vec::new();
    for dy in -r_px..=r_px {
        for dx in -r_px..=r_px {
            if (dx * dx + dy * dy) as f64 <= r2 {
                offs.push((dx, dy));
            }
        }
    }
    let w = shape.width;
    let h = shape.height;
    let mut out = vec![0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            if !valid[i] {
                continue;
            }
            let mut sum = 0f64;
            let mut cnt = 0u32;
            for &(dx, dy) in &offs {
                let nx = x as i64 + dx;
                let ny = y as i64 + dy;
                if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                    continue;
                }
                let j = ny as usize * w + nx as usize;
                if valid[j] {
                    sum += raw[j] as f64;
                    cnt += 1;
                }
            }
            out[i] = if cnt > 0 {
                (sum / cnt as f64) as f32
            } else {
                raw[i]
            };
        }
    }
    out
}
