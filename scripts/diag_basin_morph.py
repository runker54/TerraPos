# -*- coding: utf-8 -*-
"""诊断方案三坝子形态缺陷(固定路径): 对象级 长宽比/内切宽度/包围度/内部起伏"""
import numpy as np
from osgeo import gdal
from scipy import ndimage

DEM = r"G:\tif_features\county_feature\hhgq\dem.tif"
BASIN = r"E:\zcode_worker\Topographic\rust\target\hhgq_out\basin_mask.tif"

# 直接抽稀读取到 ~25m(避免 537MB 全量加载崩溃)
ds_d = gdal.Open(DEM); W, H = ds_d.RasterXSize, ds_d.RasterYSize
step = 5
print("读 DEM(抽稀)...", flush=True)
dem25 = ds_d.GetRasterBand(1).ReadAsArray(0, 0, W, H, W // step, H // step)
print("读坝子(抽稀)...", flush=True)
ds_b = gdal.Open(BASIN)
basin25 = ds_b.GetRasterBand(1).ReadAsArray(0, 0, W, H, W // step, H // step) == 1
res = 25.0

print("标记连通域...", flush=True)
lab, n = ndimage.label(basin25, structure=np.ones((3, 3), bool))
print(f"盆底对象 {n} 个(25m 层)")
sizes = ndimage.sum(basin25, lab, range(1, n + 1))

print("逐对象统计...", flush=True)
rows = []
for k in range(1, n + 1):
    m = lab == k
    a_cells = int(m.sum())
    if a_cells < 40:  # <2.5 公顷跳过细节
        continue
    # 内切宽度: 逐步腐蚀测最大内切圆半径
    md = ndimage.distance_transform_edt(m)
    w_in = 2 * md.max() * res
    # 长宽比: 面积 / (内切圆直径)^2 的倒数意义弱, 用外接盒长宽比
    ys, xs = np.where(m)
    box_ratio = (ys.max() - ys.min() + 1) / max(1, (xs.max() - xs.min() + 1))
    # 包围度: 膨胀环带 P75 - 对象中位数
    ring = ndimage.binary_dilation(m, iterations=8) & ~m  # 200m
    if ring.sum() == 0:
        continue
    ring_p75 = np.percentile(dem25[ring], 75)
    obj_med = np.median(dem25[m])
    surround = ring_p75 - obj_med
    # 内部起伏 P95-P5
    inner = np.percentile(dem25[m], 95) - np.percentile(dem25[m], 5)
    rows.append((k, a_cells * 625 / 1e6, w_in, box_ratio, surround, inner))

rows.sort(key=lambda r: -r[1])
print(f"{'编号':>4} {'面积km²':>8} {'内切宽m':>8} {'盒长宽比':>8} {'环带高差m':>9} {'内部起伏m':>8}")
for r in rows[:25]:
    print(f"{r[0]:>4} {r[1]:>8.2f} {r[2]:>8.0f} {r[3]:>8.2f} {r[4]:>9.1f} {r[5]:>8.1f}")

arr = np.array(rows)
if len(arr):
    print("\n汇总(面积>=2.5公顷对象):")
    print(f"  数量 {len(arr)}, 总面积 {arr[:,1].sum():.1f} km²")
    print(f"  内切宽度: 中位 {np.median(arr[:,2]):.0f}m, <250m 的占 {100*(arr[:,2]<250).mean():.0f}%")
    print(f"  环带高差(包围度): 中位 {np.median(arr[:,4]):.1f}m, <10m 的占 {100*(arr[:,4]<10).mean():.0f}%")
    print(f"  内部起伏: 中位 {np.median(arr[:,5]):.1f}m, >15m 的占 {100*(arr[:,5]>15).mean():.0f}%")
    long_narrow = (arr[:,3] > 5) | (arr[:,2] < 250)
    print(f"  疑似窄长谷带(盒比>5 或 内切宽<250m): {long_narrow.sum()} 个, 面积 {arr[long_narrow,1].sum():.1f} km² ({100*arr[long_narrow,1].sum()/arr[:,1].sum():.0f}%)")
