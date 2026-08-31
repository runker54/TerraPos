# -*- coding: utf-8 -*-
"""最终对比: geomorphon 坡位 vs legacy 坡位(分块读)"""
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

MINE = r"E:\zcode_worker\Topographic\rust\target\hhgq_out\slope_position.tif"
LEGACY = r"G:\tif_features\county_feature\hhgq\slopeposition.tif"

def read_all(path):
    ds = gdal.Open(path)
    W, H = ds.RasterXSize, ds.RasterYSize
    bh = (H + 7) // 8
    rows = []
    for b in range(8):
        yoff = b * bh
        ys = min(bh, H - yoff)
        rows.append(ds.GetRasterBand(1).ReadAsArray(0, yoff, W, ys))
    ds = None
    return np.vstack(rows)

a1 = read_all(MINE)
a2 = read_all(LEGACY)
agree = 100.0 * (a1 == a2).sum() / a1.size
for name, a in [("geomorphon", a1), ("legacy", a2)]:
    v, c = np.unique(a, return_counts=True)
    print(name, {int(k): round(100.0 * cc / a.size, 1) for k, cc in zip(v, c)})
print(f"逐像元一致率: {agree:.1f}%")
