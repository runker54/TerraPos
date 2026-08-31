# -*- coding: utf-8 -*-
"""
用 arcpy 在 25m DEM 上做水文分析（喀斯特区：只填浅洼，保留深洼/落水洞）
输出: work/acc25_arcpy.tif (汇流累积量)
河网阈值判断由 Python 端完成
"""
import os
import arcpy
from arcpy import env
from arcpy.sa import Fill, FlowDirection, FlowAccumulation

WORK = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'work')
DEM25 = os.path.join(WORK, 'dem25.tif')
env.overwriteOutput = True
env.workspace = WORK
arcpy.CheckOutExtension('Spatial')

print('Fill(15m)...', flush=True)
fill_r = Fill(DEM25, 15)
print('FlowDirection...', flush=True)
fdir_r = FlowDirection(fill_r, 'NORMAL')
print('FlowAccumulation...', flush=True)
acc_r = FlowAccumulation(fdir_r, None, 'FLOAT')
acc_path = os.path.join(WORK, 'acc25_arcpy.tif')
acc_r.save(acc_path)
print('已保存:', acc_path, flush=True)
arcpy.CheckInExtension('Spatial')
