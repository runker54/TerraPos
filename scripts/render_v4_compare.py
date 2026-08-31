# -*- coding: utf-8 -*-
"""v3(TPI) vs v4(OBIA) 坝子对比图(固定路径)"""
import numpy as np
from osgeo import gdal
from PIL import Image, ImageDraw, ImageFont

BASE = r"E:\zcode_worker\Topographic"
ds = gdal.Open(r"G:\tif_features\county_feature\hhgq\dem.tif")
W, H = ds.RasterXSize, ds.RasterYSize
step = 4
w4, h4 = W // step, H // step
dem = ds.GetRasterBand(1).ReadAsArray(0, 0, W, H, w4, h4)
res = ds.GetGeoTransform()[1] * step
ds = None
print("DEM 抽稀", dem.shape, flush=True)
gy, gx = np.gradient(dem, res)
slope = np.arctan(np.hypot(gx, gy))
aspect = np.arctan2(-gx, -gy)
az, alt = np.radians(315), np.radians(45)
hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
hs8 = (np.clip(hs, 0, 1) * 255).astype(np.uint8)
base = np.dstack([hs8, hs8, hs8])
del dem, gy, gx, slope, aspect, hs

ds = gdal.Open(BASE + r"\rust\target\hhgq_out\basin_mask.tif")
b3 = ds.GetRasterBand(1).ReadAsArray(0, 0, W, H, w4, h4) == 1
ds = None
ds = gdal.Open(BASE + r"\rust\target\basin_obia\basin_obia.tif")
b4 = ds.GetRasterBand(1).ReadAsArray(0, 0, W, H, w4, h4) == 1
ds = None
print("坝子读取完成", flush=True)

p3 = base.copy()
p3[b3] = (52, 178, 229)
p4 = base.copy()
p4[b4] = (52, 222, 120)
del b3, b4, base
sep = np.full((8, w4, 3), 30, np.uint8)
img = Image.fromarray(np.vstack([p3, sep, p4]))
del p3, p4
d = ImageDraw.Draw(img)
f = ImageFont.truetype("C:/Windows/Fonts/msyh.ttc", 26)
d.text((14, 10), "v3 TPI负地形: 65.4km2/150个 (24%不被山环绕, 20%内部起伏>15m 混入)", fill=(255, 255, 80), font=f)
d.text((14, h4 + 20), "v4 OBIA对象级: 76.1km2/116个 (面积+平坦+包围+宽度 四重检验)", fill=(120, 255, 150), font=f)
img.thumbnail((1500, 3000))
img.save(BASE + r"\docs\basin_v4_compare.png")
print("saved", img.size)
