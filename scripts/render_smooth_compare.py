# -*- coding: utf-8 -*-
"""平滑效果对比: 我方 SDF 版 vs Manba 参数版(shapely buffer) 局部放大"""
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from PIL import Image, ImageDraw, ImageFont

DEM = r"G:\tif_features\county_feature\hhgq\dem.tif"
HIS_SHP = r"G:\soil_shp_data\hhgq\river\alluvial_areas_final_500_5_5_5_5000_50.shp"
MINE = r"E:\zcode_worker\Topographic\rust\target\hhgq_out\basin_mask.tif"

ds = rasterio.open(DEM)
W, H = ds.width, ds.height
g0 = ds.transform
# 局部窗口: 取最大坝子质心附近 6x4 km
mine_full = rasterio.open(MINE).read(1)
from scipy import ndimage
lab, n = ndimage.label(mine_full == 1)
if n:
    sizes = ndimage.sum(mine_full == 1, lab, range(1, n + 1))
    k = int(np.argmax(sizes)) + 1
    ys, xs = np.where(lab == k)
    cy, cx = int(ys.mean()), int(xs.mean())
else:
    cy, cx = H // 2, W // 2
del mine_full, lab
half = 3000  # 6km 窗(5m)
row0, row1 = max(0, cy - half), min(H, cy + half)
col0, col1 = max(0, cx - half), min(W, cx + half)
rw, rh = col1 - col0, row1 - row0
win_transform = g0 * Affine.translation(col0, row0)
print(f"局部窗 {rw}x{rh} @({col0},{row0})", flush=True)

dem_win = ds.read(1, window=rasterio.windows.Window(col0, row0, rw, rh)).astype(np.float32)
res = g0.a
gy, gx = np.gradient(dem_win, res)
slope = np.arctan(np.hypot(gx, gy))
aspect = np.arctan2(-gx, -gy)
az, alt = np.radians(315), np.radians(45)
hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
hs8 = (np.clip(hs, 0, 1) * 255).astype(np.uint8)
del dem_win, gy, gx, slope, aspect, hs

gdf = gpd.read_file(HIS_SHP)
his = rasterize([(g, 1) for g in gdf.geometry], out_shape=(rh, rw), transform=win_transform, fill=0, dtype="uint8") == 1
del gdf
mine = rasterio.open(MINE).read(1, window=rasterio.windows.Window(col0, row0, rw, rh)) == 1

p1 = np.dstack([hs8, hs8, hs8]); p1[mine] = (52, 178, 229)
p2 = np.dstack([hs8, hs8, hs8]); p2[his] = (255, 120, 60)
sep = np.full((6, rw, 3), 30, np.uint8)
img = Image.fromarray(np.vstack([p1, sep, p2]))
d = ImageDraw.Draw(img)
f = ImageFont.truetype("C:/Windows/Fonts/msyh.ttc", 30)
d.text((16, 12), "我方(SDF 平滑重构)", fill=(120, 220, 255), font=f)
d.text((16, rh + 18), "原始(arcpy 矢量 buffer 平滑)", fill=(255, 200, 120), font=f)
img.thumbnail((1400, 1400))
img.save(r"E:\zcode_worker\Topographic\docs\smooth_compare.png")
print("saved docs/smooth_compare.png", img.size)
