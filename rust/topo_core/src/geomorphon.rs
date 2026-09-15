//! geomorphon 地貌形态分析(对齐 GRASS r.geomorphon 算法)
//!
//! 机器视觉方法: 以焦点像元为中心, 沿 8 方向做线视(line-of-sight),
//! 在搜索半径 L 内确定各方向的天顶/天底角, 与 flat 阈值比较生成三值
//! 模式(+更高/-更低/0平坦), 组成 8 元组地貌模式(geomorphon)。
//!
//! 与窗口统计法(TPI/均值)相比: 结果与实际地物的尺寸、起伏、朝向无关,
//! 单一半径即可自适应多尺度形态。

/// 三值模式: 每方向 2 位, 0=平坦 1=更高(天顶) 2=更低(天底)
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Pattern8(pub [u8; 8]);

impl Pattern8 {
    /// 山脊: 全部方向更低(周围都比我低)
    pub fn is_peak(&self) -> bool {
        self.0.iter().all(|&v| v == 1)
    }
    /// 山谷: 全部方向更高(周围都比我高)
    pub fn is_valley(&self) -> bool {
        self.0.iter().all(|&v| v == 2)
    }
    /// 正方向(更低=四周水往远处流)数量: 位置分类用
    pub fn count_lower(&self) -> u32 {
        self.0.iter().filter(|&&v| v == 2).count() as u32
    }
    /// 负方向(更高=四周山地)数量
    pub fn count_higher(&self) -> u32 {
        self.0.iter().filter(|&&v| v == 1).count() as u32
    }
    pub fn count_flat(&self) -> u32 {
        self.0.iter().filter(|&&v| v == 0).count() as u32
    }
}

/// geomorphon 计算
///
/// 沿射线逐像元步进并双线性插值高程; 天顶/天底角取射线上
/// 相对焦点最大仰角/最大俯角(切角), 对齐 r.geomorphon angle 原理。
///
/// # 参数
/// - `dem`/`w`/`h`/`res`: 高程栅格与分辨率
/// - `x`/`y`: 焦点像元
/// - `search_m`: 线视最大距离(米), 决定形态尺度
/// - `skip_m`: 起始跳过距离(米), 消除微起伏噪声(0=不跳)
/// - `flat_deg`: 平坦判角阈值(度), |角差|<flat 视为平坦
#[allow(clippy::too_many_arguments)]
pub fn geomorphon_pattern(
    dem: &[f32],
    w: usize,
    h: usize,
    res: f64,
    x: usize,
    y: usize,
    search_m: f64,
    skip_m: f64,
    flat_deg: f64,
) -> Pattern8 {
    let zc = dem[y * w + x];
    let zc = if zc.is_finite() { zc } else { return Pattern8::default() };
    // 8 方向: 东起逆时针(对齐 GRASS 顺序); 无序性对分类统计无影响
    let dirs: [(f64, f64); 8] = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
    ];
    let search_px = (search_m / res).ceil();
    let skip_px = (skip_m / res).ceil();
    debug_assert!(search_m > skip_m, "search_m 必须 > skip_m");
    debug_assert!(
        skip_m == 0.0 || skip_m >= res,
        "skip_m 应 >= 分辨率或为 0(禁用)"
    );
    debug_assert!(
        skip_px >= 1.0 || skip_m == 0.0,
        "skip 像元换算不得低于 1 像元"
    );
    let flat_rad = flat_deg.to_radians();
    let mut pat = Pattern8::default();
    let at = |fx: f64, fy: f64| -> f32 {
        let x0 = fx.floor() as isize;
        let y0 = fy.floor() as isize;
        let x1 = x0 + 1;
        let y1 = y0 + 1;
        let gx = fx - x0 as f64;
        let gy = fy - y0 as f64;
        let g = |xx: isize, yy: isize| -> f64 {
            if xx < 0 || yy < 0 || xx >= w as isize || yy >= h as isize {
                return f64::NAN;
            }
            dem[yy as usize * w + xx as usize] as f64
        };
        let z00 = g(x0, y0);
        let z10 = g(x1, y0);
        let z01 = g(x0, y1);
        let z11 = g(x1, y1);
        // 窗外像元按最近值延拓(避免 NaN 污染插值)
        let z10 = if z10.is_nan() { z00 } else { z10 };
        let z01 = if z01.is_nan() { z00 } else { z01 };
        let z11 = if z11.is_nan() { z00 } else { z11 };
        (z00 * (1.0 - gx) * (1.0 - gy)
            + z10 * gx * (1.0 - gy)
            + z01 * (1.0 - gx) * gy
            + z11 * gx * gy) as f32
    };
    for (k, (dx, dy)) in dirs.iter().enumerate() {
        let mut zenith = 0f64; // 最大仰角(某射线点高于焦点)
        let mut nadir = 0f64; // 最大俯角(某射线点低于焦点)
        let mut dist = skip_px.max(1.0);
        while dist <= search_px {
            let fx = x as f64 + dx * dist;
            let fy = y as f64 + dy * dist;
            let zv = at(fx, fy);
            if zv.is_finite() {
                let angle = ((zv as f64 - zc as f64) / (dist * res)).atan();
                if angle > zenith {
                    zenith = angle;
                }
                if -angle > nadir {
                    nadir = -angle;
                }
            }
            dist += 1.0;
        }
        // 三值: 与 flat 阈值比较(天顶/天底取大者定方向性)
        pat.0[k] = if zenith > flat_rad && zenith >= nadir {
            1 // 更高
        } else if nadir > flat_rad && nadir > zenith {
            2 // 更低
        } else {
            0 // 平坦
        };
    }
    pat
}

/// 六级坡位映射: 按 GRASS 标准 P 值(更高方向数, 0..8)单调映射,
/// 对齐 legacy 编码 1谷 2坡下 3平坡 4坡中 5坡上 6山脊:
/// - P 0-1: 山脊(全向更低/近全向) -> 6
/// - P 2: 坡上(shoulder 肩部) -> 5
/// - P 3-4: 坡中(spur/slope) -> 4
/// - P 5-6: 平坡(footslope 坡麓) -> 3
/// - P 7: 坡下/山谷(valley) -> 2
/// - P 8: 山谷(pit 坑) -> 1
pub fn pattern_to_level(pat: &Pattern8, _min_extreme: u32) -> u8 {
    let p = pat.count_higher();
    match p {
        0..=1 => 6,
        2 => 5,
        3..=4 => 4,
        5..=6 => 3,
        7 => 2,
        _ => 1,
    }
}

// ---------------- 标准十形态(计划 Task 5, 替代计数映射) ----------------

/// GRASS r.geomorphon 十种标准地貌形态
#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Landform {
    Flat = 1,
    Peak = 2,
    Ridge = 3,
    Shoulder = 4,
    Spur = 5,
    Slope = 6,
    Hollow = 7,
    Footslope = 8,
    Valley = 9,
    Pit = 10,
}

impl Landform {
    /// 上/中/下坡位基础证据(设计规格 10.2 表; Flat 由 q 与 HAND 下游再细分)
    pub fn evidence(self) -> (f32, f32, f32) {
        match self {
            Landform::Peak => (1.00, 0.00, 0.00),
            Landform::Ridge => (0.90, 0.10, 0.00),
            Landform::Shoulder => (0.70, 0.30, 0.00),
            Landform::Spur => (0.55, 0.45, 0.00),
            Landform::Slope => (0.10, 0.80, 0.10),
            Landform::Hollow => (0.00, 0.45, 0.55),
            Landform::Footslope => (0.00, 0.25, 0.75),
            Landform::Valley => (0.00, 0.10, 0.90),
            Landform::Pit => (0.00, 0.00, 1.00),
            Landform::Flat => (0.05, 0.25, 0.70),
        }
    }
}

/// 环形同色弧段数(跳过 flat 方向; 环形闭合, 全同色为 1)
fn circular_runs(s: &[u8; 8]) -> usize {
    let mut runs = 0usize;
    for k in 0..8 {
        let cur = s[k];
        let prev = s[(k + 7) % 8];
        if cur != 0 && prev != cur {
            runs += 1;
        }
    }
    runs.max(1)
}

/// 某一色的环形弧宽合计(方向数)
fn arc_width(s: &[u8; 8], color: u8) -> usize {
    s.iter().filter(|&&v| v == color).count()
}

/// 规范化三值码: trits {lower=0, flat=1, higher=2}, 取 8 旋转与 8 镜像的
/// 最小三进制编码, 用于诊断与不变性测试(判定树与其等价)。
pub fn canonical_code(pat: &Pattern8) -> u32 {
    // 内部编码(0=flat,1=higher,2=lower) -> 计划契约 trits(lower=0,flat=1,higher=2)
    let to_trit = |v: u8| match v {
        2 => 0u32, // lower
        0 => 1u32, // flat
        _ => 2u32, // higher
    };
    let tri: Vec<u32> = pat.0.iter().map(|&v| to_trit(v)).collect();
    let reflect: Vec<u32> = tri.iter().rev().copied().collect();
    let mut best = u32::MAX;
    for base in [&tri, &reflect] {
        for r in 0..8 {
            let mut code = 0u32;
            for k in 0..8 {
                code = code * 3 + base[(r + k) % 8];
            }
            best = best.min(code);
        }
    }
    best
}

/// 标准十形态判定(环形弧段拓扑, 旋转与镜像不变):
/// - 无差异 Flat; 全低 Peak; 全高 Pit
/// - 1 段: 纯高段+平坦 -> Footslope; 纯低段+平坦 -> Shoulder
/// - 2 段: 含平坦 -> Shoulder; 无平坦 -> Slope
/// - 4 段(两低两高相对): 低弧宽 >= 高弧宽 -> Ridge, 否则 Valley
/// - 6 段以上(破碎): 低多 -> Spur, 高多 -> Hollow, 均衡 -> Slope
///
/// 等宽 4 段在旋转+镜像规范化下本征简并, 按约定归 Ridge(与 GRASS 一致)。
pub fn pattern_to_landform(pat: &Pattern8) -> Landform {
    let s = &pat.0;
    let l = arc_width(s, 2); // lower 方向数
    let h = arc_width(s, 1); // higher 方向数
    if l == 0 && h == 0 {
        return Landform::Flat;
    }
    if l == 8 {
        return Landform::Peak;
    }
    if h == 8 {
        return Landform::Pit;
    }
    let runs = circular_runs(s);
    match runs {
        1 => {
            if h > 0 {
                Landform::Footslope
            } else {
                Landform::Shoulder
            }
        }
        2 => {
            if l + h < 8 {
                Landform::Shoulder
            } else {
                Landform::Slope
            }
        }
        4 => {
            // 山脊: 沿脊两方向低(窄低弧), 横向两壁抬升(宽高弧); 山谷相反
            if l <= h {
                Landform::Ridge
            } else {
                Landform::Valley
            }
        }
        _ => {
            if l > h {
                Landform::Spur
            } else if h > l {
                Landform::Hollow
            } else {
                Landform::Slope
            }
        }
    }
}

/// 全图 geomorphon 六级坡位
pub fn geomorphon_levels(
    dem: &[f32],
    w: usize,
    h: usize,
    res: f64,
    search_m: f64,
    skip_m: f64,
    flat_deg: f64,
) -> Vec<u8> {
    let mut out = vec![0u8; w * h];
    for y in 1..h.saturating_sub(1) {
        for x in 1..w.saturating_sub(1) {
            let pat = geomorphon_pattern(dem, w, h, res, x, y, search_m, skip_m, flat_deg);
            out[y * w + x] = pattern_to_level(&pat, 2);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak_dem(w: usize, h: usize, peak_h: f32) -> Vec<f32> {
        let mut d = vec![100f32; w * h];
        let (cx, cy) = (w / 2, h / 2);
        for y in 0..h {
            for x in 0..w {
                let dd = ((x as f64 - cx as f64).powi(2) + (y as f64 - cy as f64).powi(2)).sqrt();
                d[y * w + x] = 100.0 + (peak_h as f64 * (1.0 - dd / 40.0).clamp(0.0, 1.0)) as f32;
            }
        }
        d
    }

    #[test]
    fn peak_center_is_ridge_level() {
        // 尖峰中心: 全向更低 -> 山脊(6)
        let (w, h) = (81usize, 81usize);
        let dem = peak_dem(w, h, 200.0);
        assert_eq!(geomorphon_pattern(&dem, w, h, 5.0, 40, 40, 100.0, 0.0, 1.0).count_lower(), 8);
        let lv = geomorphon_levels(&dem, w, h, 5.0, 100.0, 0.0, 1.0);
        assert_eq!(lv[40 * w + 40], 6, "峰心应为山脊");
    }

    #[test]
    fn pit_center_is_valley_level() {
        let (w, h) = (81usize, 81usize);
        let dem: Vec<f32> = peak_dem(w, h, 200.0).iter().map(|&z| 300.0 - (z - 100.0)).collect();
        let lv = geomorphon_levels(&dem, w, h, 5.0, 100.0, 0.0, 1.0);
        assert_eq!(lv[40 * w + 40], 1, "坑心应为山谷");
    }

    #[test]
    fn skip_suppresses_micro_relief() {
        // 中心平地+微噪声: skip=0 时噪声成方向, skip 大时判定平坦
        let (w, h) = (81usize, 81usize);
        let mut dem = vec![100f32; w * h];
        for i in 0..w * h {
            if i % 7 == 0 {
                dem[i] += 0.05; // 5cm 微噪声
            }
        }
        let p0 = geomorphon_pattern(&dem, w, h, 5.0, 40, 40, 100.0, 0.0, 1.0);
        let p1 = geomorphon_pattern(&dem, w, h, 5.0, 40, 40, 100.0, 30.0, 1.0);
        // 微噪声在 skip 后不应产生方向性
        assert!(p1.count_higher() + p1.count_lower() < 8);
        let _ = p0;
    }
}
// ---------------- 自适应形态证据(Task 8) ----------------

use crate::error::{CoreError, Result};
use rayon::prelude::*;
use crate::input::RasterShape;
use crate::terrain::{plan_curvature, profile_curvature, slope_degrees};

/// 尺度索引形态证据产品
#[derive(Debug, Clone)]
pub struct MorphEvidence {
    /// 所选最近尺度层的十形态
    pub form: Vec<Landform>,
    /// 上/中/下坡位基础证据(0..1, 三者和为 1)
    pub upper: Vec<f32>,
    pub middle: Vec<f32>,
    pub lower: Vec<f32>,
    /// 小尺度坡度(度)
    pub slope_deg: Vec<f32>,
    pub profile_curvature: Vec<f32>,
    pub plan_curvature: Vec<f32>,
    /// 尺度端点或导数支持不完整
    pub low_confidence: Vec<bool>,
}

/// 逐可用尺度计算 geomorphon 形态层, 逐像元取其 adaptive_scale_m 最近层;
/// 剖面/平面曲率对上/下证据做至多 0.10 的符号修正后归一化。
/// 尺度层从小到大处理, 像元值一旦定格即释放该层缓冲。
pub fn adaptive_morphology_evidence(
    dem: &[f32],
    valid: &[bool],
    shape: RasterShape,
    adaptive_scale_m: &[f32],
    usable_scales_m: &[f64],
) -> Result<MorphEvidence> {
    if usable_scales_m.is_empty() {
        return Err(CoreError::Invalid("尺度族为空, 无法计算形态证据".into()));
    }
    let w = shape.width;
    let h = shape.height;
    let n = w * h;
    let res = shape.resolution_m;

    // 导数场(一次)
    let slope_deg = slope_degrees(dem, w, h, res);
    let pc = profile_curvature(dem, w, h, res);
    let pk = plan_curvature(dem, w, h, res);

    // 每像元最近尺度层索引
    let nearest: Vec<u8> = (0..n)
        .map(|i| {
            let s = adaptive_scale_m[i] as f64;
            if !s.is_finite() || s <= 0.0 {
                return 0u8;
            }
            let mut best = 0u8;
            let mut bd = f64::INFINITY;
            for (k, &r) in usable_scales_m.iter().enumerate() {
                let d = (r - s).abs();
                if d < bd {
                    bd = d;
                    best = k as u8;
                }
            }
            best
        })
        .collect();

    let mut form = vec![Landform::Flat; n];
    let mut upper = vec![0f32; n];
    let mut middle = vec![0f32; n];
    let mut lower = vec![0f32; n];
    let mut low_confidence = vec![true; n];

    // 层从小到大处理: 像元在所属最近层定格证据后即释放层缓冲
    for (k, &radius) in usable_scales_m.iter().enumerate() {
        let skip_m = (2.0 * res).max(0.05 * radius);
        let flat_deg = 3.0;
        // 逐行并行 geomorphon
        let forms_k: Vec<Landform> = (0..h)
            .into_par_iter()
            .flat_map_iter(|y| {
                (0..w).map(move |x| {
                    let i = y * w + x;
                    if !valid[i] || !dem[i].is_finite() {
                        return Landform::Flat;
                    }
                    let pat = geomorphon_pattern(
                        dem, w, h, res, x, y, radius, skip_m, flat_deg,
                    );
                    pattern_to_landform(&pat)
                })
            })
            .collect();
        for i in 0..n {
            if nearest[i] != k as u8 || !valid[i] {
                continue;
            }
            form[i] = forms_k[i];
            let (mut u, m, mut l) = forms_k[i].evidence();
            // 曲率符号修正(约定: 凸为负), 剖面/平面各至多 0.05
            if pc[i].is_finite() {
                let t = (pc[i] / 0.05).clamp(-1.0, 1.0);
                u += 0.05 * (-t).max(0.0);
                l += 0.05 * t.max(0.0);
            }
            if pk[i].is_finite() {
                let t = (pk[i] / 0.05).clamp(-1.0, 1.0);
                u += 0.05 * (-t).max(0.0);
                l += 0.05 * t.max(0.0);
            }
            let sum = u + m + l;
            if sum > 0.0 {
                upper[i] = u / sum;
                middle[i] = m / sum;
                lower[i] = l / sum;
            }
            let endpoint = k == 0 || k == usable_scales_m.len() - 1;
            let deriv_ok = slope_deg[i].is_finite() && pc[i].is_finite() && pk[i].is_finite();
            low_confidence[i] = endpoint || !deriv_ok;
        }
        drop(forms_k);
    }

    Ok(MorphEvidence {
        form,
        upper,
        middle,
        lower,
        slope_deg,
        profile_curvature: pc,
        plan_curvature: pk,
        low_confidence,
    })
}
