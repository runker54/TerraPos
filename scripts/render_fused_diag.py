# -*- coding: utf-8 -*-
"""融合版(fused_5) vs Manba 参数版 叠加诊断图"""
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from PIL import Image, ImageDraw, ImageFont

DEM = r"G:\tif_features\county_feature\hhgq\dem.tif"
HIS_SHP = r"G:\soil_shp_data\hhgq\river\alluvial_areas_final_500_5_5_5_5000_50.shp"
MINE = r"E:\zcode_worker\Topographic\rust\target\basin_fused\fused_5.tif"
RIVER_SHP = r"G:\soil_shp_data\hhgq\hhgq_sd_polygon.shp"

ds = rasterio.open(DEM)
W, H = ds.width, ds.height
g0 = ds.transform
step = 4
w, h = W // step, H // step
transform = Affine(g0.a * step, g0.b, g0.c, g0.d, g0.e * step, g0.f)
dem = ds.read(1, out_shape=(h, w)).astype(np.float32)

gdf = gpd.read_file(HIS_SHP)
his = rasterize([(g, 1) for g in gdf.geometry], out_shape=(h, w), transform=transform, fill=0, dtype="uint8") == 1
del gdf
mine = rasterio.open(MINE).read(1, out_shape=(h, w)) == 1

# 地类河流(1101/1102/1103/1107)栅格化
sd = gpd.read_file(RIVER_SHP)
riv = sd[sd["DLBM"].isin(["1101", "1102", "1103", "1107"])]
river = rasterize([(g, 1) for g in riv.geometry], out_shape=(h, w), transform=transform, fill=0, dtype="uint8") == 1
del sd, riv
print(f"地类河流要素: {river.sum()} 像元({20}m)", flush=True)

res = g0.a * step
gy, gx = np.gradient(dem, res)
slope = np.arctan(np.hypot(gx, gy))
aspect = np.arctan2(-gx, -gy)
az, alt = np.radians(315), np.radians(45)
hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
hs8 = (np.clip(hs, 0, 1) * 255).astype(np.uint8)
rgb = np.dstack([hs8, hs8, hs8])
rgb[his & mine] = (60, 200, 90)     # 绿: 双方
rgb[mine & ~his] = (70, 150, 255)   # 蓝: 仅我
rgb[his & ~mine] = (255, 90, 70)    # 红: 仅他
# 河流描亮
rgb[river] = np.where(rgb[river] == 0, (240, 240, 200), rgb[river]).astype(np.uint8) if False else rgb[river]
rgb[river] = np.clip(rgb[river].astype(np.int32) + 60, 0, 255).astype(np.uint8)

img = Image.fromarray(rgb)
d = ImageDraw.Draw(img)
f = ImageFont.truetype("C:/Windows/Fonts/msyh.ttc", 22)
d.text((12, 10), "绿=双方 蓝=仅我方 红=仅原始 亮色=地类河流", fill=(255, 255, 60), font=f)
img.thumbnail((1500, 1500))
img.save(r"E:\zcode_worker\Topographic\docs\basin_fused_diag.png")
print("saved docs/basin_fused_diag.png")
