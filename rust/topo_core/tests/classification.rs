//! 分类契约测试: 模糊隶属度、置信度与确定性坡位(Task 11 后扩展组合与清理)。

// ---------------- Task 9: 模糊上中下隶属度 ----------------

use topo_core::geomorphon::{Landform, MorphEvidence};
use topo_core::slope_position::{classify_slope_positions, SlopePosition};

fn neutral_morph(n: usize) -> MorphEvidence {
    MorphEvidence {
        form: vec![Landform::Flat; n],
        upper: vec![0.0; n],
        middle: vec![0.0; n],
        lower: vec![0.0; n],
        slope_deg: vec![5.0; n],
        profile_curvature: vec![0.0; n],
        plan_curvature: vec![0.0; n],
        low_confidence: vec![false; n],
    }
}

/// q=0.1/0.5/0.9 处隶属度排序正确; 中心平局取中部
#[test]
fn fuzzy_membership_orders_by_relative_position() {
    let qs: Vec<f32> = vec![0.1, 0.5, 0.9];
    let n = qs.len();
    let morph = neutral_morph(n);
    let m = classify_slope_positions(&qs, &morph, &vec![true; n]).unwrap();
    assert_eq!(m.raw[0], SlopePosition::Lower, "q=0.1 应为下部");
    assert_eq!(m.raw[1], SlopePosition::Middle, "q=0.5 中心应为中部(平局取中)");
    assert_eq!(m.raw[2], SlopePosition::Upper, "q=0.9 应为上部");
}

/// 0.38/0.62 过渡带两侧隶属度连续(无跳变)
#[test]
fn fuzzy_membership_is_continuous_at_transitions() {
    let qs: Vec<f32> = vec![0.37, 0.39, 0.61, 0.63];
    let morph = neutral_morph(qs.len());
    let m = classify_slope_positions(&qs, &morph, &vec![true; qs.len()]).unwrap();
    for k in [0usize, 2] {
        let d_upper = (m.upper[k] - m.upper[k + 1]).abs();
        let d_lower = (m.lower[k] - m.lower[k + 1]).abs();
        assert!(d_upper < 0.15, "上隶属度跨 0.62 跳变: {d_upper}");
        assert!(d_lower < 0.15, "下隶属度跨 0.38 跳变: {d_lower}");
    }
}

/// 形态证据可扭转接近的平局: 中部位置 + 强峰证据 -> 上部
#[test]
fn morphology_resolves_near_tie() {
    let n = 1;
    let mut morph = neutral_morph(n);
    morph.upper[0] = 5.0; // 强峰证据
    let m = classify_slope_positions(&[0.5f32], &morph, &[true]).unwrap();
    assert_eq!(m.raw[0], SlopePosition::Upper, "强峰证据应扭转中部主导");
}

/// 置信度 = 最大隶属度 - 次大隶属度; 低置信折扣 0.6 不改变排序
#[test]
fn confidence_and_low_confidence_discount() {
    let n = 2;
    let qs = [0.9f32, 0.9];
    let mut morph = neutral_morph(n);
    let m = classify_slope_positions(&qs, &morph, &vec![true; n]).unwrap();
    let mut s: Vec<f32> = vec![m.upper[0], m.middle[0], m.lower[0]];
    s.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let expect = s[0] - s[1];
    assert!((m.confidence[0] - expect).abs() < 1e-4, "置信度错误");
    // 低置信(端点尺度) -> 乘 0.6
    morph.low_confidence[1] = true;
    let m2 = classify_slope_positions(&qs, &morph, &vec![true; n]).unwrap();
    assert!((m2.confidence[1] - m.confidence[1] * 0.6).abs() < 1e-4);
    assert_eq!(m2.raw[1], m.raw[0], "折扣不得改变排序");
}

/// q 无效(锚点缺失)时仅凭形态证据分类且置信度打折
#[test]
fn missing_geometry_falls_back_to_morphology() {
    let n = 1;
    let mut morph = neutral_morph(n);
    morph.lower[0] = 1.0;
    morph.form[0] = Landform::Pit;
    let m = classify_slope_positions(&[f32::NAN], &morph, &[true]).unwrap();
    assert_eq!(m.raw[0], SlopePosition::Lower, "仅形态证据时应为下部");
    assert!(m.confidence[0] <= 0.6, "几何缺失置信度应打折");
}
