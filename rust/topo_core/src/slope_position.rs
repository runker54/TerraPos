//! 模糊上/中/下坡位: 连续隶属度、置信度与确定性坡位。
//!
//! 主位置量 q(融合水平/垂直相对位置)经 Logistic/高斯隶属函数生成
//! 基础隶属度, 形态证据做线性修正; 无 q 时退化为纯形态证据分类。

use crate::error::Result;
use crate::geomorphon::MorphEvidence;

/// 三级坡位(内部序; 组合层映射到最终编码 3-8)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlopePosition {
    Lower = 1,
    Middle = 2,
    Upper = 3,
}

/// 隶属度产品
#[derive(Debug, Clone)]
pub struct Memberships {
    pub upper: Vec<f32>,
    pub middle: Vec<f32>,
    pub lower: Vec<f32>,
    /// 确定性坡位(最大隶属度)
    pub raw: Vec<SlopePosition>,
    /// 置信度 = 最大隶属度 - 次大隶属度(低置信 *0.6)
    pub confidence: Vec<f32>,
}

#[inline]
fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}

/// 分类模糊坡位。q 无效(锚点缺失)时基础隶属度取均匀 1/3, 仅凭形态证据。
pub fn classify_slope_positions(
    relative_position: &[f32],
    morph: &MorphEvidence,
    valid: &[bool],
) -> Result<Memberships> {
    let n = relative_position.len();
    let mut upper = vec![0f32; n];
    let mut middle = vec![0f32; n];
    let mut lower = vec![0f32; n];
    let mut raw = vec![SlopePosition::Middle; n];
    let mut confidence = vec![0f32; n];
    for i in 0..n {
        if !valid[i] {
            continue;
        }
        let q = relative_position[i];
        let (u0, m0, l0) = if q.is_finite() {
            (
                sigmoid((q - 0.62) / 0.08),
                (-((q - 0.50) / 0.22).powi(2)).exp(),
                sigmoid((0.38 - q) / 0.08),
            )
        } else {
            (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
        };
        let u = u0 + 0.20 * morph.upper[i];
        let m = m0 + 0.15 * morph.middle[i];
        let l = l0 + 0.20 * morph.lower[i];
        let sum = u + m + l;
        let (u, m, l) = if sum > 0.0 { (u / sum, m / sum, l / sum) } else { (u, m, l) };
        // 确定性: 最大者; 精确平局优先级 中部 > 下部 > 上部
        let (r, top, second) = if u > m && u > l {
            (SlopePosition::Upper, u, m.max(l))
        } else if m > l {
            (SlopePosition::Middle, m, if u > l { u } else { l })
        } else {
            (SlopePosition::Lower, l, m.max(u))
        };
        upper[i] = u;
        middle[i] = m;
        lower[i] = l;
        raw[i] = r;
        let conf = top - second;
        confidence[i] = if morph.low_confidence[i] || !q.is_finite() {
            conf * 0.6
        } else {
            conf
        };
    }
    Ok(Memberships {
        upper,
        middle,
        lower,
        raw,
        confidence,
    })
}
