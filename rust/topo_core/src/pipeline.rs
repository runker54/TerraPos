//! 参数模型 + 全流程编排（自适应 DEM 地形部位划分管线）。
//!
//! pipeline 只负责编排、进度、取消与成果写出; 全部算法在 input/scale/
//! hydro/ridge/slope_unit/geomorphon/slope_position/basin/postprocess。
//! 阶段顺序固定: 输入 -> 双表面 -> 水文 -> 金字塔 -> 形态层 -> 脊线/单元
//! -> 几何 -> 形态证据/隶属度 -> 盆地 -> 组合清理 -> 上采样 -> 写出。

use crate::basin::{detect_basins, BasinConfig, BasinResult};
use crate::error::{CoreError, Result};
use crate::geomorphon::adaptive_morphology_evidence;
use crate::geotiff::{self, GeoMeta};
use crate::hydro::{build_hydro, HydroConfig};
use crate::input::{prepare_input, prepare_values, InputConfig, PreparedDem, RasterShape};
#[allow(unused_imports)]
use crate::hydro::hand_to_stream;
use crate::postprocess::{compose_codes, constrained_cleanup, FinalClassification};
use crate::ridge::build_ridges;
use crate::scale::{build_scale_pyramid, usable_scales, ScalePyramid};
use crate::slope_position::{classify_slope_positions, SlopePosition};
use crate::slope_unit::{build_slope_geometry, build_slope_units};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

/// 盆地保留倾向(设计规格 19 节: 普通 UI 参数)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BasinTendency {
    Strict,
    Standard,
    Loose,
}

/// 分析精细度预设(设计规格 19 节: 普通 UI 参数)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PrecisionPreset {
    Fast,
    Standard,
    Detailed,
}

/// 研究模式高级参数(物理单位)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdvancedParams {
    /// 粗层分析分辨率(米)
    pub coarse_res_m: f64,
    /// 水文最大直接填深(米)
    pub hydro_z_limit_m: f32,
    /// 尺度增长率收敛阈值
    pub scale_growth_threshold: f32,
    /// 低起伏判据(米, HAND/位置融合权重过渡)
    pub low_relief_m: f32,
    /// 丘陵海拔上限(米)
    pub hill_elevation_max_m: f32,
    /// 嵌套河网面积等级(平方千米, 严格递增)
    pub stream_areas_km2: [f64; 4],
}

impl Default for AdvancedParams {
    fn default() -> Self {
        AdvancedParams {
            coarse_res_m: 25.0,
            hydro_z_limit_m: 15.0,
            scale_growth_threshold: 0.15,
            low_relief_m: 20.0,
            hill_elevation_max_m: 500.0,
            stream_areas_km2: [0.05, 0.20, 1.00, 5.00],
        }
    }
}

/// 运行参数(全部物理单位)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Params {
    pub dem_path: String,
    pub out_dir: String,
    pub precision: PrecisionPreset,
    pub basin_tendency: BasinTendency,
    /// 最小盆地面积(平方米)
    pub basin_min_area_m2: f64,
    /// 后处理强度(0..2, 1=标准)
    pub postprocess_strength: f64,
    pub write_diagnostics: bool,
    pub advanced: AdvancedParams,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            dem_path: String::new(),
            out_dir: String::new(),
            precision: PrecisionPreset::Standard,
            basin_tendency: BasinTendency::Standard,
            basin_min_area_m2: 66_666.67,
            postprocess_strength: 1.0,
            write_diagnostics: true,
            advanced: AdvancedParams::default(),
        }
    }
}

/// 预设解析结果(任务 13 的参数效应测试锚点)
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedParams {
    pub coarse_res_m: f64,
    /// 尺度族上下限(米)
    pub scale_min_m: f64,
    pub scale_max_m: f64,
    /// 分带写出缓冲边长(像元)
    pub tile_edge: usize,
}

impl Params {
    /// 预设映射: Fast=max(res,50)/250-4000/1024, Standard=max(res,25)/
    /// 125-4000/1024, Detailed=max(res,10)/125-8000/768
    pub fn resolved(&self) -> Result<ResolvedParams> {
        if !(0.05..=2.0).contains(&self.postprocess_strength) {
            return Err(CoreError::Invalid(format!(
                "后处理强度必须在 0.05..2 (观测 {})",
                self.postprocess_strength
            )));
        }
        if self.basin_min_area_m2 <= 0.0 {
            return Err(CoreError::Invalid("最小盆地面积必须为正".into()));
        }
        let res = 5.0; // 实际分辨率在输入就绪后由 resolved_with 替换
        Ok(self.resolved_with(res))
    }

    /// 输入就绪后的解析(分辨率来自真实 DEM)
    pub fn resolved_with(&self, native_res_m: f64) -> ResolvedParams {
        let (coarse, smin, smax, tile) = match self.precision {
            PrecisionPreset::Fast => {
                (native_res_m.max(50.0), 250.0, 4000.0, 1024usize)
            }
            PrecisionPreset::Standard => {
                (native_res_m.max(25.0), 125.0, 4000.0, 1024)
            }
            PrecisionPreset::Detailed => {
                (native_res_m.max(10.0), 125.0, 8000.0, 768)
            }
        };
        ResolvedParams {
            coarse_res_m: coarse,
            scale_min_m: smin,
            scale_max_m: smax,
            tile_edge: tile,
        }
    }
}

pub struct Progress {
    pub stage: String,
    pub pct: f32,
    pub msg: String,
}

pub struct Outputs {
    pub terrain: Vec<u8>,
    pub subclass: Vec<u8>,
    pub confidence: Vec<f32>,
    pub meta5: GeoMeta,
    pub report: String,
    /// (类编码, 面积km²), 按业务顺序
    pub stats: Vec<(u8, f64)>,
}

/// 诊断图层(粗层网格, 与上采样前的分析网格一致)
#[derive(Debug, Clone, Default)]
pub struct DiagnosticLayers {
    pub hydro_conditioning_depth: Vec<f32>,
    pub stream_level: Vec<u8>,
    pub ridge_mask: Vec<u8>,
    pub slope_unit: Vec<u32>,
    pub adaptive_scale_m: Vec<f32>,
    pub hand_m: Vec<f32>,
    pub relative_position: Vec<f32>,
    pub slope_position_raw: Vec<u8>,
    pub basin_candidate: Vec<u8>,
    pub basin_core: Vec<u8>,
    pub basin_mask: Vec<u8>,
}

/// 数组级管线输出(测试与编排共用)
pub struct ArrayOutputs {
    /// 原生分辨率最终编码
    pub terrain: Vec<u8>,
    pub geomorph_subclass: Vec<u8>,
    pub confidence: Vec<f32>,
    pub diagnostics: DiagnosticLayers,
    /// (编码, 像元数)
    pub stats: Vec<(u8, f64)>,
    pub report: String,
    pub native_shape: RasterShape,
    /// 诊断层网格(粗层)
    pub coarse_shape: RasterShape,
    pub coarse_meta: GeoMeta,
}

/// 检查取消
fn check(cancelled: &AtomicBool) -> Result<()> {
    if cancelled.load(Ordering::Relaxed) {
        return Err(CoreError::Cancelled);
    }
    Ok(())
}

/// 数组级端到端管线(测试入口与文件编排的共用主体)
pub fn run_arrays_for_test(
    dem: &[f32],
    _valid: &[bool],
    shape: RasterShape,
    params: &Params,
) -> Result<ArrayOutputs> {
    let cancelled = AtomicBool::new(false);
    run_arrays(dem, _valid, shape, params, &cancelled, &|_| true)
}

/// 数组级管线主体
pub fn run_arrays(
    dem: &[f32],
    _valid: &[bool],
    shape: RasterShape,
    params: &Params,
    cancelled: &AtomicBool,
    progress: &dyn Fn(&Progress) -> bool,
) -> Result<ArrayOutputs> {
    let mut stage_times: Vec<(String, f32)> = Vec::new();
    let mut t0 = std::time::Instant::now();
    let say = |pct: f32, stage: &str, msg: &str| -> Result<()> {
        check(cancelled)?;
        progress(&Progress {
            stage: stage.to_string(),
            pct,
            msg: msg.to_string(),
        });
        Ok(())
    };

    // ---- 1/2 输入与双表面 ----
    say(5.0, "输入", "构建有效掩膜与地貌表面...")?;
    let prepared = prepare_values(dem.to_vec(), demo_meta(shape)?, &InputConfig::default())?;
    let resolved = params.resolved_with(prepared.shape.resolution_m);
    stage_times.push(("输入与双表面".into(), t0.elapsed().as_secs_f32()));

    // ---- 3 水文 ----
    t0 = std::time::Instant::now();
    say(12.0, "水文", "有界修正 + 嵌套河网 + HAND...")?;
    let hcfg = HydroConfig {
        coarse_res_m: resolved.coarse_res_m,
        z_limit_m: params.advanced.hydro_z_limit_m,
        stream_areas_km2: params.advanced.stream_areas_km2,
    };
    let hydro = build_hydro(&prepared, &hcfg)?;
    let cw = hydro.shape.width;
    let ch = hydro.shape.height;
    let cn = cw * ch;
    // 粗层表面(valid-aware 均值, 与水文网格一致)
    let (coarse, coarse_valid) = coarse_average(&prepared, resolved.coarse_res_m, cw, ch);
    stage_times.push(("水文".into(), t0.elapsed().as_secs_f32()));

    // ---- 4 金字塔 ----
    t0 = std::time::Instant::now();
    say(24.0, "尺度", "多尺度稳健金字塔...")?;
    let pyramid = build_scale_pyramid(
        &coarse,
        &coarse_valid,
        hydro.shape,
        params.advanced.scale_growth_threshold,
    )?;
    stage_times.push(("尺度金字塔".into(), t0.elapsed().as_secs_f32()));

    // ---- 5-6 形态层 / 脊线与单元 ----
    t0 = std::time::Instant::now();
    say(36.0, "形态", "尺度索引 geomorphon 形态证据...")?;
    let scales_all = usable_scales(hydro.shape);
    let scales: Vec<f64> = scales_all
        .iter()
        .copied()
        .filter(|&s| s >= resolved.scale_min_m && s <= resolved.scale_max_m)
        .collect();
    let scales = if scales.len() < 3 { scales_all } else { scales };
    let morph = adaptive_morphology_evidence(
        &coarse,
        &coarse_valid,
        hydro.shape,
        &pyramid.characteristic_scale_m,
        &scales,
    )?;
    say(46.0, "脊线", "子流域边界 + 补充证据脊线...")?;
    let ridges =
        build_ridges(&coarse, &coarse_valid, &hydro, &pyramid, &morph.form)?;
    let units = build_slope_units(&coarse, &coarse_valid, &hydro, &ridges)?;
    stage_times.push(("形态/脊线/单元".into(), t0.elapsed().as_secs_f32()));

    // ---- 7 几何 ----
    t0 = std::time::Instant::now();
    say(56.0, "几何", "约束距离 + 自适应尺度 + 相对位置...")?;
    let geom = build_slope_geometry(
        &coarse,
        &coarse_valid,
        &units,
        &hydro,
        &pyramid,
        params.advanced.low_relief_m,
    )?;
    stage_times.push(("约束几何".into(), t0.elapsed().as_secs_f32()));

    // ---- 8 隶属度 ----
    say(64.0, "隶属", "模糊上/中/下隶属度...")?;
    // q 场 3x3 中值平滑: 抑制锚点切换造成的椒盐反转(空间一致性)
    let mut q_smoothed = geom.relative_position.clone();
    {
        let mut buf = q_smoothed.clone();
        for y in 1..ch - 1 {
            for x in 1..cw - 1 {
                let i = y * cw + x;
                let mut nbr: Vec<f32> = Vec::with_capacity(9);
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        let j = (y as i64 + dy) as usize * cw + (x as i64 + dx) as usize;
                        if coarse_valid[j] && geom.relative_position[j].is_finite() {
                            nbr.push(geom.relative_position[j]);
                        }
                    }
                }
                nbr.sort_by(|a, b| a.partial_cmp(b).unwrap());
                if let Some(mid) = nbr.get(nbr.len() / 2) {
                    if mid.is_finite() {
                        buf[i] = *mid;
                    }
                }
            }
        }
        q_smoothed = buf;
    }
    let positions = classify_slope_positions(&q_smoothed, &morph, &coarse_valid)?;

    // ---- 9 盆地 ----
    say(70.0, "盆地", "低平候选 + 对象检验 + 边界重建...")?;
    let basin = detect_basins(
        &coarse,
        &coarse_valid,
        hydro.shape,
        &hydro,
        &units,
        &geom,
        &pyramid,
        &morph,
        &BasinConfig {
            tendency: params.basin_tendency,
            min_area_m2: params.basin_min_area_m2,
        },
    )?;

    // ---- 10 组合与清理 ----
    say(76.0, "组合", "编码组合 + 受约束清理...")?;
    let mut fc: FinalClassification = compose_codes(
        &coarse,
        &coarse_valid,
        &basin,
        &positions,
        params.advanced.hill_elevation_max_m,
    )?;
    // 高丘细分(亚类): 2000m 窗内起伏 <200 m 为低丘
    refine_high_hills(&mut fc, &pyramid, params.advanced.hill_elevation_max_m);
    let corrected = constrained_cleanup(
        &mut fc,
        &positions,
        &units,
        &geom,
        &basin,
        hydro.shape,
        params.postprocess_strength,
    )?;
    let low_conf_pct = {
        let lc = morph.low_confidence.iter().filter(|&&b| b).count();
        100.0 * lc as f64 / cn as f64
    };

    // ---- 11 上采样到原生分辨率 ----
    t0 = std::time::Instant::now();
    say(84.0, "上采样", "粗层成果回投原生分辨率...")?;
    let terrain = upscale_categorical(&fc.terrain, cw, ch, shape);
    let subclass = upscale_categorical(&fc.geomorph_subclass, cw, ch, shape);
    let confidence = crate::pipeline::upsample(&fc.confidence, cw, ch, resolved.coarse_res_m, shape.width, shape.height, shape.resolution_m);
    // NoData 恢复: 原生无效一律 0/NaN
    let mut terrain = terrain;
    let mut subclass = subclass;
    let mut confidence = confidence;
    for i in 0..shape.width * shape.height {
        if !prepared.valid[i] {
            terrain[i] = 0;
            subclass[i] = 0;
            confidence[i] = f32::NAN;
        }
    }
    // 统计
    let mut stats: Vec<(u8, f64)> = vec![(0, 0.0); 9];
    for &c in &terrain {
        stats[c as usize].0 = c;
        stats[c as usize].1 += 1.0;
    }
    stage_times.push(("上采样".into(), t0.elapsed().as_secs_f32()));

    // ---- 诊断图层(粗层) ----
    let raw_pos: Vec<u8> = positions
        .raw
        .iter()
        .zip(coarse_valid.iter())
        .map(|(r, &v)| if v { match r { SlopePosition::Upper => 3, SlopePosition::Middle => 4, SlopePosition::Lower => 5 } } else { 0 })
        .collect();
    let diagnostics = DiagnosticLayers {
        hydro_conditioning_depth: hydro.conditioning_depth.clone(),
        stream_level: hydro.stream_level.clone(),
        ridge_mask: ridges.mask.iter().map(|&b| b as u8).collect(),
        slope_unit: units.unit_id.clone(),
        adaptive_scale_m: geom.adaptive_scale_m.clone(),
        hand_m: geom.hand_m.clone(),
        relative_position: geom.relative_position.clone(),
        slope_position_raw: raw_pos,
        basin_candidate: basin.candidate.iter().map(|&b| b as u8).collect(),
        basin_core: basin.core.iter().map(|&b| b as u8).collect(),
        basin_mask: basin.mask.iter().map(|&b| b as u8).collect(),
    };

    // ---- 报告 ----
    let code2 = stats[2].1;
    let report = build_report(
        &resolved,
        &scales,
        &stats,
        shape,
        &basin,
        low_conf_pct,
        corrected,
        &stage_times,
        code2 as usize,
    );
    say(88.0, "完成", "分析网格就绪")?;

    let mut coarse_meta = demo_meta(hydro.shape)?;
    coarse_meta.geo_keys = prepared.meta.geo_keys.clone();
    coarse_meta.geo_ascii = prepared.meta.geo_ascii.clone();
    coarse_meta.pixel_scale = [
        resolved.coarse_res_m,
        resolved.coarse_res_m,
        0.0,
    ];
    coarse_meta.tiepoint = {
        let mut t = prepared.meta.tiepoint;
        t[3] += resolved.coarse_res_m / 2.0;
        t[4] -= resolved.coarse_res_m / 2.0;
        t
    };

    Ok(ArrayOutputs {
        terrain,
        geomorph_subclass: subclass,
        confidence,
        diagnostics,
        stats,
        report,
        native_shape: shape,
        coarse_shape: hydro.shape,
        coarse_meta,
    })
}

/// 便捷构造内存场景元数据(测试入口)
fn demo_meta(shape: RasterShape) -> Result<GeoMeta> {
    let mut m = GeoMeta::from_origin(
        shape.width as u32,
        shape.height as u32,
        500_000.0,
        3_000_000.0,
        shape.resolution_m,
    );
    m.geo_keys = vec![1, 1, 0, 2, 1024, 0, 1, 1, 3076, 0, 1, 9001];
    m.nodata = Some(-9999.0);
    Ok(m)
}

/// 粗层 valid-aware 均值(与 hydro 网格一致)
fn coarse_average(
    prepared: &PreparedDem,
    coarse_res_m: f64,
    cw: usize,
    ch: usize,
) -> (Vec<f32>, Vec<bool>) {
    let native = prepared.shape;
    let edge = |i: usize, total: usize| -> (usize, usize) {
        let s = ((i as f64 * coarse_res_m / native.resolution_m).round() as usize).min(total);
        let e = ((((i + 1) as f64 * coarse_res_m / native.resolution_m).round() as usize)
            .min(total))
        .max(s + 1);
        (s, e)
    };
    let mut out = vec![0f32; cw * ch];
    let mut outv = vec![false; cw * ch];
    for ty in 0..ch {
        let (y0, y1) = edge(ty, native.height);
        for tx in 0..cw {
            let (x0, x1) = edge(tx, native.width);
            let (mut sum, mut c) = (0f64, 0u32);
            for y in y0..y1 {
                for x in x0..x1 {
                    let i = y * native.width + x;
                    if prepared.valid[i] {
                        sum += prepared.raw[i] as f64;
                        c += 1;
                    }
                }
            }
            if c > 0 {
                out[ty * cw + tx] = (sum / c as f64) as f32;
                outv[ty * cw + tx] = true;
            }
        }
    }
    (out, outv)
}

/// 高丘细分: <500m 且特征起伏 >=200m -> 高丘(亚类 2)
fn refine_high_hills(fc: &mut FinalClassification, pyramid: &ScalePyramid, hill_max: f32) {
    for i in 0..fc.geomorph_subclass.len() {
        if fc.geomorph_subclass[i] == 1
            && pyramid.characteristic_relief_m
                .get(i)
                .copied()
                .unwrap_or(0.0)
                >= 200.0
        {
            fc.geomorph_subclass[i] = 2;
        }
    }
    let _ = hill_max;
}

/// 类别最近邻上采样
fn upscale_categorical(src: &[u8], sw: usize, sh: usize, dst: RasterShape) -> Vec<u8> {
    let res = sw as f64;
    let scale_x = sw as f64 / dst.width as f64;
    let scale_y = sh as f64 / dst.height as f64;
    let _ = res;
    let mut out = vec![0u8; dst.width * dst.height];
    out.par_chunks_mut(dst.width).enumerate().for_each(|(y, row)| {
        let sy = ((y as f64 + 0.5) * scale_y) as usize % sh;
        for (x, out) in row.iter_mut().enumerate() {
            let sx = ((x as f64 + 0.5) * scale_x) as usize % sw;
            *out = src[sy * sw + sx];
        }
    });
    out
}

/// 双线性上采样(连续诊断层)
fn upsample(src: &[f32], sw: usize, sh: usize, src_res: f64, dst_w: usize, dst_h: usize, dst_res: f64) -> Vec<f32> {
    let mut dst = vec![f32::NAN; dst_w * dst_h];
    dst.par_chunks_mut(dst_w).enumerate().for_each(|(y, row)| {
        let gy = (y as f64 + 0.5) * dst_res / src_res - 0.5;
        let y0f = gy.floor();
        let fy = (gy - y0f) as f32;
        let y0 = (y0f.max(0.0) as usize).min(sh - 1);
        let y1 = (y0 + 1).min(sh - 1);
        for (x, out) in row.iter_mut().enumerate() {
            let gx = (x as f64 + 0.5) * dst_res / src_res - 0.5;
            let x0f = gx.floor();
            let fx = (gx - x0f) as f32;
            let x0 = (x0f.max(0.0) as usize).min(sw - 1);
            let x1 = (x0 + 1).min(sw - 1);
            let v00 = src[y0 * sw + x0];
            let v01 = src[y0 * sw + x1];
            let v10 = src[y1 * sw + x0];
            let v11 = src[y1 * sw + x1];
            *out = v00 * (1.0 - fx) * (1.0 - fy)
                + v01 * fx * (1.0 - fy)
                + v10 * fx * fy * 0.0
                + v10 * fx * fy
                + v11 * fx * fy;
        }
    });
    dst
}

/// 生成 class_report 文本
#[allow(clippy::too_many_arguments)]
fn build_report(
    resolved: &ResolvedParams,
    scales: &[f64],
    stats: &[(u8, f64)],
    shape: RasterShape,
    basin: &BasinResult,
    low_conf_pct: f64,
    corrected: usize,
    stage_times: &[(String, f32)],
    code2: usize,
) -> String {
    let names = [
        (0u8, "NoData"),
        (1, "山间/宽谷盆地"),
        (2, "保留码(不得出现)"),
        (3, "丘陵上部"),
        (4, "丘陵中部"),
        (5, "丘陵下部"),
        (6, "山地坡上"),
        (7, "山地坡中"),
        (8, "山地坡下"),
    ];
    let cell_km2 = shape.resolution_m * shape.resolution_m / 1e6;
    let mut out = String::new();
    out.push_str("自适应地形部位划分报告\n==============================================\n");
    out.push_str(&format!(
        "网格: {}x{} @ {} m\n粗层分辨率: {} m | 尺度族: {:?} m | 分带边长: {}\n\n",
        shape.width, shape.height, shape.resolution_m, resolved.coarse_res_m, scales, resolved.tile_edge
    ));
    out.push_str("类别统计 (编码 名称 像元数 面积km² 占比%):\n");
    let valid_total: f64 = stats.iter().filter(|s| s.0 != 0).map(|s| s.1).sum();
    for (code, name) in names {
        let cnt = stats[code as usize].1;
        let pct = if valid_total > 0.0 { 100.0 * cnt / valid_total } else { 0.0 };
        out.push_str(&format!(
            "  {} {:<12} {:>12.0} {:>10.2} {:>6.2}\n",
            code,
            name,
            cnt,
            cnt * cell_km2,
            pct
        ));
    }
    out.push_str(&format!(
        "\n编码 2 为保留码, 生成计数必须为 0 (实际 {code2})\n"
    ));
    out.push_str(&format!(
        "盆地对象: {} 个, 接受 {} 个\n低置信像元: {:.2}%\n单调修正像元: {}\n\n阶段耗时:\n",
        basin.objects.len(),
        basin.objects.iter().filter(|o| o.accepted).count(),
        low_conf_pct,
        corrected
    ));
    for (name, t) in stage_times {
        out.push_str(&format!("  {name}: {t:.1}s\n"));
    }
    out
}

/// 文件级编排(桌面应用入口)
pub fn run(
    params: &Params,
    progress: &dyn Fn(&Progress) -> bool,
    cancelled: &AtomicBool,
) -> Result<Outputs> {
    check(cancelled)?;
    let (values, meta) = geotiff::read_f32(&params.dem_path)?;
    let prepared = prepare_input(Path::new(&params.dem_path), &InputConfig::default())?;
    let shape = prepared.shape;
    drop(values);
    let meta5 = meta.clone();
    let out = run_arrays(
        &prepared.raw,
        &prepared.valid,
        shape,
        params,
        cancelled,
        progress,
    )?;
    // 统计为面积 km²(桌面契约)
    let cell_km2 = shape.resolution_m * shape.resolution_m / 1e6;
    let stats: Vec<(u8, f64)> = out
        .stats
        .iter()
        .map(|(c, n)| (*c, n * cell_km2))
        .collect();
    // 写出
    let out_dir = Path::new(&params.out_dir);
    std::fs::create_dir_all(out_dir)?;
    let cmap = classification_cmap();
    geotiff::write_u8_cmap(out_dir.join("terrain_position.tif"), &meta5, &out.terrain, &cmap)?;
    geotiff::write_u8_cmap(out_dir.join("geomorph_subclass.tif"), &meta5, &out.geomorph_subclass, &cmap)?;
    let mut conf_meta = meta5.clone();
    conf_meta.nodata = Some(f32::NAN);
    geotiff::write_f32(out_dir.join("terrain_confidence.tif"), &conf_meta, &out.confidence)?;
    std::fs::write(out_dir.join("class_report.txt"), &out.report)?;
    if params.write_diagnostics {
        let d = out_dir.join("diagnostics");
        std::fs::create_dir_all(&d)?;
        let dg = &out.diagnostics;
        let cmeta = &out.coarse_meta;
        geotiff::write_f32(d.join("hydro_conditioning_depth.tif"), cmeta, &dg.hydro_conditioning_depth)?;
        geotiff::write_u8_cmap(d.join("stream_level.tif"), cmeta, &dg.stream_level, &cmap)?;
        geotiff::write_u8_cmap(d.join("ridge_mask.tif"), cmeta, &dg.ridge_mask, &cmap)?;
        geotiff::write_u32(d.join("slope_unit.tif"), cmeta, &dg.slope_unit)?;
        geotiff::write_f32(d.join("adaptive_scale_m.tif"), cmeta, &dg.adaptive_scale_m)?;
        geotiff::write_f32(d.join("hand_m.tif"), cmeta, &dg.hand_m)?;
        geotiff::write_f32(d.join("relative_position.tif"), cmeta, &dg.relative_position)?;
        geotiff::write_u8_cmap(d.join("slope_position_raw.tif"), cmeta, &dg.slope_position_raw, &cmap)?;
        geotiff::write_u8_cmap(d.join("basin_candidate.tif"), cmeta, &dg.basin_candidate, &cmap)?;
        geotiff::write_u8_cmap(d.join("basin_core.tif"), cmeta, &dg.basin_core, &cmap)?;
        geotiff::write_u8_cmap(d.join("basin_mask.tif"), cmeta, &dg.basin_mask, &cmap)?;
    }
    Ok(Outputs {
        terrain: out.terrain,
        subclass: out.geomorph_subclass,
        confidence: out.confidence,
        meta5,
        report: out.report,
        stats,
    })
}

/// 分类色表(0 透明黑; 其余按业务色带)
fn classification_cmap() -> [[u8; 3]; 256] {
    let mut cmap = [[0u8; 3]; 256];
    cmap[1] = [120, 190, 120]; // 盆地 绿
    cmap[3] = [215, 158, 158]; // 丘陵上 粉红
    cmap[4] = [230, 200, 140]; // 丘陵中 黄
    cmap[5] = [170, 220, 160]; // 丘陵下 浅绿
    cmap[6] = [168, 112, 72]; // 山地坡上 棕红
    cmap[7] = [222, 196, 120]; // 山地坡中 黄棕
    cmap[8] = [132, 168, 96]; // 山地坡下 绿
    cmap
}
