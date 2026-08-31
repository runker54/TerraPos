# -*- coding: utf-8 -*-
"""我方坝子 vs Manba 原始坝区结果 对比(25m 层直接栅格化, 避免全尺寸内存)"""
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from scipy import ndimage
from PIL import Image, ImageDraw, ImageFont

DEM = r"G:\tif_features\county_feature\hhgq\dem.tif"
MINE = r"E:\zcode_worker\Topographic\rust\target\hhgq_out\basin_mask.tif"
HIS_SHP = r"G:\soil_shp_data\hhgq\river\alluvial_areas_final_result.shp"

import os
os.environ.setdefault("GDAL_CACHEMAX", "512")

w, h = 3707, 2264  # 20m 层(14831/4, 9058/4)
gt = (267000.0, 20.0, 0.0, 2866000.0, 0.0, -20.0)  # 由读取得到, 下行修正

ds = rasterio.open(DEM)
W, H = ds.width, ds.height
g0 = ds.transform
step = 4
w, h = W // step, H // step
transform = Affine(g0.a * step, g0.b, g0.c, g0.d, g0.e * step, g0.f)
dem25 = ds.read(1, out_shape=(h, w)).astype(np.float32)
ds = None
print("网格", w, "x", h, flush=True)

gdf = gpd.read_file(HIS_SHP)
print("原始坝区矢量:", len(gdf), "个面要素", flush=True)
his = rasterize(
    [(geom, 1) for geom in gdf.geometry],
    out_shape=(h, w), transform=transform, fill=0, dtype="uint8"
) == 1
del gdf
print("栅格化完成", flush=True)

ds2 = rasterio.open(MINE)
mine = ds2.read(1, out_shape=(h, w)) == 1
ds2 = None

cell = 20.0 * 20.0 / 1e6
a_his = his.sum() * cell
a_mine = mine.sum() * cell
inter = (his & mine).sum() * cell
union = (his | mine).sum() * cell
iou = inter / union if union else 0
print(f"\nManba 原始: {a_his:.2f} km2")
print(f"我方 v4.1 : {a_mine:.2f} km2")
print(f"交集 {inter:.2f} | IoU {100*iou:.1f}%")
print(f"他有我无(漏检): {(his & ~mine).sum()*cell:.2f} km2 ({100*(his&~mine).sum()/max(1,his.sum()):.0f}%)")
print(f"我有他无(多检): {(mine & ~his).sum()*cell:.2f} km2", flush=True)

st8 = np.ones((3, 3), bool)
_, n_his = ndimage.label(his, structure=st8)
_, n_mine = ndimage.label(mine, structure=st8)
sz_his = ndimage.sum(his, _, range(1, n_his + 1))
sz_mine = ndimage.sum(mine, _, range(1, n_mine + 1))
print(f"\n原始对象 {n_his} 个, 面积中位 {np.median(sz_his)*cell*100:.1f} 公顷, 最大 {sz_his.max()*cell*100:.0f} 公顷")
print(f"我方对象 {n_mine} 个, 面积中位 {np.median(sz_mine)*cell*100:.1f} 公顷, 最大 {sz_mine.max()*cell*100:.0f} 公顷")

res = 20.0
gy, gx = np.gradient(dem25, res)
slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, -gy)
az, alt = np.radians(315), np.radians(45)
hs = np.sin(alt)*np.cos(slope) + np.cos(alt)*np.sin(slope)*np.cos(az-aspect)
hs8 = (np.clip(hs, 0, 1)*255).astype(np.uint8)
rgb = np.dstack([hs8, hs8, hs8])
rgb[his & mine] = (60, 200, 90)
rgb[mine & ~his] = (70, 150, 255)
rgb[his & ~mine] = (255, 90, 70)
img = Image.fromarray(rgb)
d = ImageDraw.Draw(img)
f = ImageFont.truetype("C:/Windows/Fonts/msyh.ttc", 22)
d.text((12, 10), "绿=双方一致  蓝=仅我方  红=仅原始(Manba)", fill=(255, 255, 60), font=f)
img.thumbnail((1500, 1500))
img.save(r"E:\zcode_worker\Topographic\docs\basin_vs_legacy.png")
print("saved docs/basin_vs_legacy.png")
