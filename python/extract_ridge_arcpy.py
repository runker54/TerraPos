# -*- coding: utf-8 -*-
"""
用 arcpy 在 25m DEM 上做负地形水文分析（标准山脊线提取法）
负地形 = 区域最高程 - DEM；负地形的汇流高值线即正地形的山脊/分水线网络
输出: work/acc_neg25.tif (负地形汇流累积量)
"""
import os
import arcpy
from arcpy import env
from arcpy.sa import Fill, FlowDirection, FlowAccumulation, Raster

WORK = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'work')
DEM25 = os.path.join(WORK, 'dem25.tif')
env.overwriteOutput = True
env.workspace = WORK
arcpy.CheckOutExtension('Spatial')

print('构建负地形...', flush=True)
import numpy as np
_arr = arcpy.RasterToNumPyArray(DEM25, nodata_to_value=np.nan)
zmax = float(np.nanmax(_arr))
del _arr
neg_r = Raster(DEM25)
neg_r = zmax - neg_r
print(f'区域最高程 {zmax}m', flush=True)

print('负地形全填洼(连通山脊线)...', flush=True)
fill_r = Fill(neg_r)
print('负地形流向...', flush=True)
fdir_r = FlowDirection(fill_r, 'NORMAL')
print('负地形汇流累积...', flush=True)
acc_r = FlowAccumulation(fdir_r, None, 'FLOAT')
acc_path = os.path.join(WORK, 'acc_neg25.tif')
acc_r.save(acc_path)
print('已保存:', acc_path, flush=True)
