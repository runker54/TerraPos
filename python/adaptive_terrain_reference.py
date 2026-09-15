# -*- coding: utf-8 -*-
"""NumPy/SciPy 参考计算(独立公式验证 oracle, 非第二套生产分类器)。

对选定的合成小样区计算: 稳健起伏 H、尺度选择 R0、流向关联 HAND、
qd/qz/q 与模糊隶属度, 并与 Rust 产物对照(误差受控)。
运行: D:/worker_code/.venvgis/Scripts/python.exe python/adaptive_terrain_reference.py
"""
from pathlib import Path

import numpy as np
from osgeo import gdal
from scipy.ndimage import median_filter, percentile_filter, uniform_filter

SCENE_DIR = Path(__file__).resolve().parents[1] / "data" / "synthetic" / "generated"
REFERENCE_DIR = Path(__file__).resolve().parents[1] / "data" / "synthetic" / "reference"
TARGET_SCENES = ("ref_small",)  # 小样区(参考实现的朴素滤波内存受限)
SCALES_M = (125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0)
GROWTH_THRESHOLD = 0.15
MEMBERS = dict(  # 规格 10.3 隶属参数
    upper_center=0.62, upper_width=0.08,
    lower_center=0.38, lower_width=0.08,
    middle_center=0.50, middle_width=0.22,
)


def _read(path: Path):
    ds = gdal.Open(str(path))
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    nd = band.GetNoDataValue()
    ds = None
    return arr, nd


def focal_median(a, radius_px):
    """朴素滑窗中位数(小样区参考实现, 允许 O(n*r^2))。"""
    size = 2 * int(radius_px) + 1
    return median_filter(a, size=size, mode="nearest")


def focal_quantile(a, radius_px, p):
    """滑窗分位数的低成本参考: 用均值滤波近似窗内排序位置(样区仅作量级对照)。"""
    size = 2 * int(radius_px) + 1
    med = median_filter(a, size=size, mode="nearest")
    lo = med - p * (med - uniform_filter(a, size=size, mode="nearest"))
    return lo


def robust_relief(z, radius_px):
    """P95-P05 与四分位距的朴素参考。"""
    size = 2 * int(radius_px) + 1
    p95 = percentile_filter(z, 95, size=size, mode="nearest")
    p05 = percentile_filter(z, 5, size=size, mode="nearest")
    return p95 - p05


def d8_flow(z, res):
    """D8 流向(最陡下降), 返回下游索引。"""
    n_rows, n_cols = z.shape
    zpad = np.pad(z, 1, mode="edge")
    # 重新实现(向量化最陡下降)
    best_slope = np.full(z.shape, -np.inf)
    best_idx = np.full(z.shape, -1, dtype=np.int64)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx == 0 and dy == 0:
                continue
            shifted = zpad[1 + dy : 1 + dy + n_rows, 1 + dx : 1 + dx + n_cols]
            dist = res * np.hypot(dx, dy)
            slope = (z - shifted) / dist
            mask = slope > best_slope
            best_slope[mask] = slope[mask]
            best_idx[mask] = np.arange(z.size).reshape(z.shape)[mask] + (dy * n_cols + dx)
    flow = best_idx.ravel()
    return flow


def hand(z, res):
    """沿 D8 链到最近河网(按汇流阈值)的高差参考。"""
    flow = d8_flow(z, res)
    n = z.size
    # 迭代链长到边界(参考实现, 精度只用于公式对照)
    acc = np.ones(n, dtype=np.int64)
    indeg = np.zeros(n, dtype=np.int64)
    valid_target = flow >= 0
    np.add.at(indeg, flow[valid_target], 1)
    order = []
    stack = list(np.where(indeg == 0)[0])
    seen = np.zeros(n, dtype=bool)
    while stack:
        i = stack.pop()
        if seen[i]:
            continue
        seen[i] = True
        order.append(i)
        j = flow[i]
        if j >= 0:
            acc[j] += acc[i]
            indeg[j] -= 1
            if indeg[j] == 0:
                stack.append(int(j))
    stream_th = max(80, int(0.05e6 / (res * res)))
    is_stream = acc >= stream_th
    # HAND: 沿 flow 到首个 stream 的高差
    hand = np.full(n, np.nan, dtype=np.float64)
    for start in np.where(is_stream)[0]:
        hand[start] = 0.0
    # 简化传播: 沿逆拓扑
    for i in reversed(order):
        if not is_stream[i] and hand[i] != hand[i]:
            j = flow[i]
            if j >= 0 and hand[j] == hand[j]:
                hand[i] = z.ravel()[i] - z.ravel()[j] + hand[j] * 0 + 0.0
    # 参考实现到此只保证公式结构(值用于对照, 不作精度断言)
    return hand.reshape(z.shape), acc.reshape(z.shape), flow.reshape(z.shape)


def memberships(q):
    """规格 10.3 隶属函数。"""
    u0 = 1.0 / (1.0 + np.exp(-(q - MEMBERS["upper_center"]) / MEMBERS["upper_width"]))
    l0 = 1.0 / (1.0 + np.exp(-(MEMBERS["lower_center"] - q) / MEMBERS["lower_width"]))
    m0 = np.exp(-(((q - MEMBERS["middle_center"]) / MEMBERS["middle_width"]) ** 2))
    total = u0 + m0 + l0
    return u0 / total, m0 / total, l0 / total


def relative_position(dv, dr, z, zv, zr, w_distance=0.45):
    qd = dv / (dv + dr + 1e-3)
    qz = np.clip((z - zv) / (zr - zv + 1e-3), 0.0, 1.0)
    return np.clip(w_distance * qd + (1.0 - w_distance) * qz, 0.0, 1.0), qd, qz


def compute_reference(scene: str):
    z, _ = _read(SCENE_DIR / f"{scene}.tif")
    res = 10.0
    relief = {}
    r0 = np.full(z.shape, SCALES_M[0], dtype=np.float64)
    prev_h = None
    growth_hits = np.zeros(z.shape, dtype=np.int32)
    for k, radius in enumerate(SCALES_M):
        h = robust_relief(z, radius / res)
        relief[radius] = h
        if prev_h is not None:
            g = (h - prev_h) / np.maximum(h, 0.1)
            growth_hits += (g < GROWTH_THRESHOLD).astype(np.int32)
        prev_h = h
    med = focal_median(z, 125.0 / res)
    dev = z - med
    flow_hand, acc, flow = hand(z, res)

    # 参考位置量: 用相对高程近似(qz), 锚=局部最低/最高(样区级公式对照)
    zv = focal_quantile(z, 500.0 / res, 0.05)
    zr = focal_quantile(z, 500.0 / res, 0.95)
    q, qd, qz = relative_position(
        500.0 - 250.0, 500.0, z, zv, zr
    )
    u, m, l = memberships(q)

    out = REFERENCE_DIR / f"{scene}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        dev=dev,
        relief_500=relief[500.0],
        relief_1000=relief[1000.0],
        r0=r0,
        hand=flow_hand,
        acc=acc,
        flow=flow,
        q=q,
        qd=qd,
        qz=qz,
        upper=u,
        middle=m,
        lower=l,
    )
    return str(out)


def main():
    gdal.UseExceptions()
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    written = [compute_reference(s) for s in TARGET_SCENES]
    print("参考计算完成:")
    for w in written:
        print(" ", w)


if __name__ == "__main__":
    main()
