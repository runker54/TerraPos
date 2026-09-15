//! 确定性合成 DEM 夹具与断言辅助(无随机成分)。
//!
//! 所有高程函数一律使用像元中心米坐标: `x=(col+0.5)*res`, `y=(row+0.5)*res`,
//! 保证同一场景在不同分辨率下解析一致, 供尺度稳定性测试复用。

use topo_core::geotiff::GeoMeta;
use topo_core::input::{prepare_values, RasterShape};

/// 构造带 GeoKey 目录的内存 GeoMeta(模型类型 + 线性单位可指定),
/// 供输入契约测试表达投影米制/地理度/投影英尺等元数据组合。
pub fn meta_with_keys(
    width: u32,
    height: u32,
    res: f64,
    model_type: u16,
    linear_units: u16,
) -> GeoMeta {
    let mut m = GeoMeta::from_origin(width, height, 500_000.0, 3_000_000.0, res);
    // GeoKey 目录: 头 4 SHORT(版本1, 修订1.0, 键数2) + 每键 4 SHORT(内联值)
    m.geo_keys = vec![
        1, 1, 0, 2, //
        1024, 0, 1, model_type, //
        3076, 0, 1, linear_units,
    ];
    m
}

/// 平面斜坡: 自南向北线性上升(坡降 0.08), 南缘谷底、北缘山脊。
pub fn planar_slope(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let slope = 0.08;
    let mut z = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = (row as f64 + 0.5) * res;
        for _col in 0..width {
            z.push((900.0 + slope * y) as f32);
        }
    }
    (z, RasterShape { width, height, resolution_m: res })
}

/// V 形谷: 中央谷底、两侧线性上升至边缘山脊; 谷底向北微倾保证排水。
pub fn v_valley(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let mut z = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = (row as f64 + 0.5) * res;
        for col in 0..width {
            let x = (col as f64 + 0.5) * res - cx;
            z.push((800.0 + 0.10 * x.abs() + 0.004 * y) as f32);
        }
    }
    (z, RasterShape { width, height, resolution_m: res })
}

/// 宽坦盆地: 圆形平底加环状山脊围限, 盆底沿 x 微倾提供排水出口。
pub fn broad_basin(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z: Vec<f32> = (0..height)
        .flat_map(|row| {
            (0..width).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r = x.hypot(y);
                let floor = 800.0 + 0.0005 * x;
                let rim = ((r - 700.0) / 250.0).clamp(0.0, 1.0);
                (floor + 120.0 * rim * rim) as f32
            })
        })
        .collect();
    (z, RasterShape { width, height, resolution_m: res })
}

/// 闭合洼地: 外环高地包围中央低地且无地表出口(喀斯特/山间闭合盆地原型)。
pub fn closed_pit(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z: Vec<f32> = (0..height)
        .flat_map(|row| {
            (0..width).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r = x.hypot(y);
                let v = if r < 400.0 {
                    810.0 + 2.0 * (r / 400.0) * (r / 400.0)
                } else if r < 800.0 {
                    let t = (r - 400.0) / 400.0;
                    812.0 + 60.0 * t * t
                } else {
                    872.0 + 0.08 * (r - 800.0)
                };
                v as f32
            })
        })
        .collect();
    (z, RasterShape { width, height, resolution_m: res })
}

/// 锥形山丘: 中央峰顶(850 m)加周围平原(700 m)。
pub fn conical_hill(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let cx = width as f64 * res * 0.5;
    let cy = height as f64 * res * 0.5;
    let z: Vec<f32> = (0..height)
        .flat_map(|row| {
            (0..width).map(move |col| {
                let x = (col as f64 + 0.5) * res - cx;
                let y = (row as f64 + 0.5) * res - cy;
                let r = x.hypot(y);
                let v = 700.0 + 150.0 * (1.0 - r / 600.0).max(0.0);
                v as f32
            })
        })
        .collect();
    (z, RasterShape { width, height, resolution_m: res })
}

/// 嵌套山脊: 沿 x 方向周期性脊—谷交替(谷底 850 m、脊顶 930 m, 波长 640 m),
/// 整体向北缓倾 0.8% 提供汇流方向。
pub fn nested_ridges(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    const WAVELENGTH_M: f64 = 640.0;
    let mut z = Vec::with_capacity(width * height);
    for row in 0..height {
        let y = (row as f64 + 0.5) * res;
        for col in 0..width {
            let x = (col as f64 + 0.5) * res;
            let ridge = 40.0 * (1.0 - (std::f64::consts::TAU * x / WAVELENGTH_M).cos());
            z.push((850.0 + ridge + 0.008 * y) as f32);
        }
    }
    (z, RasterShape { width, height, resolution_m: res })
}

/// 边缘 NoData: 在锥形山丘场景外圈固定 8 像元范围写入 NaN,
/// 模拟与图幅边界连通的无效区(必须保持无效, 不得内插)。
pub fn with_border_nodata(width: usize, height: usize, res: f64) -> (Vec<f32>, RasterShape) {
    let (mut z, shape) = conical_hill(width, height, res);
    let border = 8usize;
    for row in 0..height {
        for col in 0..width {
            if row < border || col < border || row + border >= height || col + border >= width {
                z[row * width + col] = f32::NAN;
            }
        }
    }
    (z, shape)
}

// ---------------- 共享测试上下文(粗层全链) ----------------

use topo_core::geomorphon::{geomorphon_pattern, pattern_to_landform, Landform};
use topo_core::hydro::{build_hydro, HydroConfig, HydroModel};
use topo_core::ridge::build_ridges;
use topo_core::slope_unit::build_slope_units;
use topo_core::scale::{build_scale_pyramid, ScalePyramid};

pub struct Context {
    pub coarse: Vec<f32>,
    pub valid: Vec<bool>,
    pub shape: RasterShape,
    pub hydro: HydroModel,
    pub pyramid: ScalePyramid,
    pub landform: Vec<Landform>,
}

fn hydro_cfg() -> HydroConfig {
    HydroConfig {
        coarse_res_m: 25.0,
        z_limit_m: 15.0,
        stream_areas_km2: [0.05, 0.20, 1.00, 5.00],
    }
}

pub fn build_context(dem: Vec<f32>, w: u32, h: u32, res: f64) -> Context {
    let prepared = prepare_values(
        dem,
        meta_with_keys(w, h, res, 1, 9001),
        &Default::default(),
    )
    .unwrap();
    let hydro = build_hydro(&prepared, &hydro_cfg()).unwrap();
    let cw = hydro.shape.width;
    let ch = hydro.shape.height;
    let coarse_res = hydro_cfg().coarse_res_m;
    let native = prepared.shape;
    // 米制边界映射, 与 hydro::downsample_bounded 完全一致
    let edge = |i: usize, total: usize| -> (usize, usize) {
        let s = ((i as f64 * coarse_res / native.resolution_m).round() as usize).min(total);
        let e = ((((i + 1) as f64 * coarse_res / native.resolution_m).round() as usize)
            .min(total))
        .max(s + 1);
        (s, e)
    };
    let mut coarse = vec![0f32; cw * ch];
    let mut valid = vec![false; cw * ch];
    for ty in 0..ch {
        let (y0, y1) = edge(ty, native.height);
        for tx in 0..cw {
            let (x0, x1) = edge(tx, native.width);
            let (mut mn, mut c) = (0f64, 0u32);
            for y in y0..y1 {
                for x in x0..x1 {
                    let i = y * native.width + x;
                    if prepared.valid[i] {
                        mn += prepared.raw[i] as f64;
                        c += 1;
                    }
                }
            }
            if c > 0 {
                coarse[ty * cw + tx] = (mn / c as f64) as f32;
                valid[ty * cw + tx] = true;
            }
        }
    }
    let pyramid =
        build_scale_pyramid(&coarse, &valid, hydro.shape, 0.15).unwrap();
    let mut landform = Vec::with_capacity(cw * ch);
    for ty in 0..ch {
        for tx in 0..cw {
            let i = ty * cw + tx;
            if !valid[i] {
                landform.push(Landform::Flat);
                continue;
            }
            let scale = pyramid.characteristic_scale_m[i].max(100.0) as f64;
            let pat = geomorphon_pattern(
                &coarse,
                cw,
                ch,
                hydro_cfg().coarse_res_m,
                tx,
                ty,
                scale,
                2.0 * hydro_cfg().coarse_res_m,
                3.0,
            );
            landform.push(pattern_to_landform(&pat));
        }
    }
    Context {
        coarse,
        valid,
        shape: hydro.shape,
        hydro,
        pyramid,
        landform,
    }
}

/// 完整链上下文: 含脊线/单元/约束几何/自适应形态证据
pub fn full_chain(
    dem: Vec<f32>,
    w: u32,
    h: u32,
    res: f64,
) -> (
    Context,
    topo_core::ridge::RidgeModel,
    topo_core::slope_unit::SlopeUnits,
    topo_core::slope_unit::SlopeGeometry,
    topo_core::geomorphon::MorphEvidence,
) {
    let ctx = build_context(dem, w, h, res);
    let ridges =
        build_ridges(&ctx.coarse, &ctx.valid, &ctx.hydro, &ctx.pyramid, &ctx.landform).unwrap();
    let units = build_slope_units(&ctx.coarse, &ctx.valid, &ctx.hydro, &ridges).unwrap();
    let geom = topo_core::slope_unit::build_slope_geometry(
        &ctx.coarse,
        &ctx.valid,
        &units,
        &ctx.hydro,
        &ctx.pyramid,
        20.0,
    )
    .unwrap();
    let scales = topo_core::scale::usable_scales(ctx.shape);
    let morph = topo_core::geomorphon::adaptive_morphology_evidence(
        &ctx.coarse,
        &ctx.valid,
        ctx.shape,
        &geom.adaptive_scale_m,
        &scales,
    )
    .unwrap();
    (ctx, ridges, units, geom, morph)
}
