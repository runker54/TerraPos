# -*- coding: utf-8 -*-
"""
西南区地形部位划分主流程 v3（HAND/RHAND 水文地形指数体系，纯DEM）

方法（第一性原理 + 国际前沿）
==========================
坡位的物理本质是"沿坡面的位置"：
  坡下 = 高于最近排水线的垂直距离小  ->  HAND  = z - z(最近河)   (Height Above Nearest Drainage,
          Rennó et al. 2008, 全球坡位/洪泛区/水文分析标准指数)
  坡上 = 低于最近分水线的垂直距离小  ->  RHAND = z(最近脊) - z   (负地形水文分析对偶, 标准脊线提取法)
  坡中 = 谷脊之间的坡腹
阈值全部为米制并按地貌亚类自适应（SUBCLASS_SLOPEPOS表），物理可解释、跨区域可迁移。
HAND/RHAND采用欧氏最近近似（25m层EDT），沟谷密集区最近河/脊即所在谷/脊，精度满足坡位分带。

山间盆地（第一性：盆地的几何本质是"开阔的低地"）：
  平缓(坡度<6°) + 显著低于周边(2km TPI<-25m) + 图斑面积>=0.5km² + 内切圆半径>=250m(剔除窄沟)
  按Manba要求宽谷盆地合并入山间盆地，编码2保留。

丘陵/山地判定（Manba口径）：
  丘陵: 海拔<500m（低丘:2km起伏<200m；高丘:200~500m）
  山地: 海拔>=500m（低山500~1000m、中山1000~3500m、高山3500~5000m、极高山>=5000m）
  本区海拔626~1714m，坡地全部为山地坡位，丘陵编码保留供其他区域。

编码（沿用calc_dxbw.py）：1山间盆地 2宽谷盆地(保留) 3丘陵上部 4丘陵中部 5丘陵下部
                          6山地坡上 7山地坡中 8山地坡下 0 NoData

依赖：
  work/dem25.tif            25m网格(gdalwarp -tr 25 25 -r average data/dem.tif work/dem25.tif)
  work/acc25_arcpy.tif      河网汇流累积(extract_river_arcpy.py)
  work/acc_neg25.tif        负地形汇流累积(extract_ridge_arcpy.py)

输出：
  output/terrain_position_5m.tif     8类地形部位栅格(带色表)
  output/geomorph_subclass_5m.tif    地貌亚类栅格(1低丘2高丘3低山4中山5高山6极高山7平坝)
  output/class_report.txt            面积统计
"""
import os
import gc
import numpy as np
import rasterio
from rasterio.warp import reproject
from rasterio.enums import Resampling as RRes
from scipy import ndimage
from osgeo import gdal

gdal.UseExceptions()

# ============================== 参数配置区 ==============================
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEM5 = os.path.join(ROOT, 'data', 'dem.tif')
ACC25 = os.path.join(ROOT, 'work', 'acc25_arcpy.tif')      # 河网汇流累积
WORK = os.path.join(ROOT, 'work')
OUT = os.path.join(ROOT, 'output')

RES_F = 5.0
RES_C = 25.0

# ---- 坡面单元分割(山体个体 watershed)与坡位判据 ----
# 以局部峰顶为种子、沿负地形分水岭分割: 每个山包/台地块一个单元, 天然"从顶到谷"完整;
# 大山坡面上的小丘包只要峰顶间距超过窗口, 即独立成单元、内部独立三等分。
# 窗口语义: WINDOW_PEAK_M = 可分辨地貌个体的最小峰顶间距(米)。
WINDOW_PEAK_M = 500.0                     # 窗口参数: 地貌个体最小间距(米)
UNIT_DEM_SMOOTH = 5                       # 峰顶提取用DEM平滑(25m像元)
UNIT_MIN_AREA = 25                        # 单元最小面积(25m像元), 小于此并入邻近单元
POS_DN = 1.0 / 3.0            # t_unit < 1/3 坡下(单元内三等分)
POS_UP = 2.0 / 3.0            # t_unit > 2/3 坡上
TERRACE_SLOPE = 2.0           # 阶地修正: 坡度<2° 且 HAND<80m 判坡下
TERRACE_HAND = 80.0

# ---- 丘陵/山地判定与亚类 ----
HILL_Z_MAX = 500.0           # 丘陵海拔上限(米)
RELIEF_SUBCLASS_WIN = 2000.0 # 丘陵亚类细分起伏度窗口(米)
RELIEF_LOW_HILL = 200.0      # 低丘相对高差上限

# ---- 水文阈值(阶地修正的HAND河网锚点) ----
ACC_RIVER_TH = 0.15e6 / (RES_C * RES_C)   # 河网集水阈值(0.15km²)

# ---- 坝区判据 ----
BASIN_SLOPE_TH = 6.0
BASIN_TPI_WIN = 2000.0
BASIN_TPI_TH = -25.0
BASIN_MIN_AREA = 500000.0    # 图斑面积下限(m²)
BASIN_CORE_HALF_W = 100.0    # 核心化半宽下限(m): 坝子全宽须>200m(Manba口径)。
                             # 实测本区候选带中位全宽86~234m, 此口径下仅真正开阔盆底可保留。
BASIN_CORE_MIN_AREA = 250000.0  # 核心化后盆底最小连片面积(m², 0.25km²=25ha)

CLASS_NAMES = {
    1: '山间盆地', 2: '宽谷盆地',
    3: '丘陵上部', 4: '丘陵中部', 5: '丘陵下部',
    6: '山地坡上', 7: '山地坡中', 8: '山地坡下',
}
CLASS_COLORS = {
    1: (51, 178, 229), 2: (102, 217, 242),
    3: (250, 217, 89), 4: (217, 237, 166), 5: (153, 199, 102),
    6: (250, 165, 60), 7: (222, 100, 50), 8: (107, 68, 35),
}
SUB_NAMES = {1: '低丘', 2: '高丘', 3: '低山', 4: '中山', 5: '高山', 6: '极高山', 7: '平坝'}
SUB_COLORS = {
    1: (180, 230, 150), 2: (110, 195, 110), 3: (250, 225, 130),
    4: (235, 170, 90), 5: (205, 110, 75), 6: (150, 65, 60), 7: (85, 185, 235),
}

# ============================== 工具函数 ==============================
def log(msg):
    print(msg, flush=True)


def focal_relief(a, wc):
    return (ndimage.maximum_filter(a, size=wc, mode='nearest')
            - ndimage.minimum_filter(a, size=wc, mode='nearest'))


def nearest_diff(z, mask):
    """z - z(最近mask像元): HAND/RHAND的欧氏近似, 返回带符号高差与距离"""
    _, (iy, ix) = ndimage.distance_transform_edt(~mask, return_indices=True)
    zv = np.where(mask, z, np.nan)
    z_near = zv[iy, ix]
    return (z - z_near).astype(np.float32)


def upsample(arr25, transform25, resampling=RRes.bilinear):
    with rasterio.open(DEM5) as ref:
        dst = np.full((ref.height, ref.width), np.nan, np.float32)
        reproject(arr25.astype(np.float32), dst,
                  src_transform=transform25, src_crs=ref.crs,
                  dst_transform=ref.transform, dst_crs=ref.crs,
                  resampling=resampling, dst_nodata=np.nan)
    return dst


def fill_nodata(a, nod):
    if nod is not None:
        a[a == nod] = np.nan
    invalid = ~np.isfinite(a)
    if invalid.any():
        _, (iy, ix) = ndimage.distance_transform_edt(invalid, return_indices=True)
        a[invalid] = a[iy[invalid], ix[invalid]]
        del iy, ix
        log(f'    nodata {100*invalid.mean():.3f}% 已填充')
    del invalid
    gc.collect()
    return a


def majority_filter_3x3(arr, classes, iters=1):
    cls_arr = np.array(classes, dtype=arr.dtype)
    lookup = np.zeros(256, np.int64)
    for i, c in enumerate(classes):
        lookup[c] = i
    for _ in range(iters):
        counts = np.stack([ndimage.uniform_filter((arr == c).astype(np.float32), size=3)
                           for c in classes])
        in_domain = np.isin(arr, cls_arr)
        idx = np.where(in_domain, lookup[arr], 0)
        orig_cnt = np.take_along_axis(counts, idx[None], axis=0)[0]
        new_arr = cls_arr[counts.argmax(axis=0)]
        keep = in_domain & (orig_cnt >= counts.max(axis=0) - 1e-6)
        new_arr[keep] = arr[keep]
        new_arr[~in_domain] = arr[~in_domain]
        arr = new_arr
        del counts
        gc.collect()
    return arr


# ============================== 主流程 ==============================
def main():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)

    for p in (DEM5, ACC25):
        assert os.path.exists(p), f'缺少 {p}'

    with rasterio.open(DEM5) as ref:
        H5, W5 = ref.height, ref.width
        transform5 = ref.transform
        crs5 = ref.crs
        dem5 = ref.read(1).astype(np.float32)
        dem5 = fill_nodata(dem5, ref.nodata)
    log(f'[0] 5m DEM 载入 {W5}x{H5}')

    # ---------- 阶段1: 25m层 HAND/RHAND 与亚类 ----------
    log('[1] 25m层 HAND/RHAND 计算...')
    with rasterio.open(os.path.join(WORK, 'dem25.tif')) as src:
        dem25 = src.read(1).astype(np.float32)
        transform25 = src.transform
        dem25 = fill_nodata(dem25, src.nodata)

    with rasterio.open(ACC25) as src:
        acc = np.where(src.read(1) == src.nodata, 0, src.read(1))
    river = acc >= ACC_RIVER_TH
    del acc
    gc.collect()
    log(f'    河网 {100*river.mean():.2f}%')

    hand25 = np.abs(nearest_diff(dem25, river))    # 高于最近河(阶地修正用)
    del river
    gc.collect()

    relief2k25 = focal_relief(dem25, int(RELIEF_SUBCLASS_WIN / RES_C))
    sub25 = np.zeros(dem25.shape, np.uint8)
    sub25[(dem25 < HILL_Z_MAX) & (relief2k25 < RELIEF_LOW_HILL)] = 1
    sub25[(dem25 < HILL_Z_MAX) & (relief2k25 >= RELIEF_LOW_HILL)] = 2
    sub25[(dem25 >= HILL_Z_MAX) & (dem25 < 1000)] = 3
    sub25[(dem25 >= 1000) & (dem25 < 3500)] = 4
    sub25[(dem25 >= 3500) & (dem25 < 5000)] = 5
    sub25[dem25 >= 5000] = 6
    present = sorted(int(s) for s in np.unique(sub25) if s > 0)
    log('    地貌亚类: ' + ', '.join(f'{s}({SUB_NAMES[s]}) {100*(sub25==s).mean():.1f}%' for s in present))

    # ---------- 坡面单元分割(峰顶种子 watershed, 窗口=峰顶最小间距) ----------
    log('[1.5] 山体单元分割(watershed)...')
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed
    dem_s = ndimage.uniform_filter(dem25, size=UNIT_DEM_SMOOTH, mode='nearest')
    wc = max(int(WINDOW_PEAK_M / RES_C), 3)
    coords = peak_local_max(dem_s, min_distance=wc, exclude_border=False)
    markers = np.zeros(dem25.shape, np.int32)
    markers[coords[:, 0], coords[:, 1]] = np.arange(1, len(coords) + 1)
    units = watershed(-dem_s, markers=markers).astype(np.int32)
    n_units = int(units.max())
    log(f'    峰顶种子 {len(coords)} 个 -> 山体单元 {n_units} 个')

    # 碎单元并入邻近单元(沿边界的最近单元): 用EDT把碎单元像元重新分配
    u_sizes = np.bincount(units.ravel())
    small_ids = np.where((u_sizes[1:] > 0) & (u_sizes[1:] < UNIT_MIN_AREA))[0] + 1
    if small_ids.size:
        small_mask = np.isin(units, small_ids)
        _, (iy_s, ix_s) = ndimage.distance_transform_edt(small_mask, return_indices=True)
        units[small_mask] = units[iy_s[small_mask], ix_s[small_mask]]
        # 被并入的单元标签可能仍引用自身(相邻全是小单元时), 再label一次紧凑化
        uniq = np.unique(units)
        remap = np.zeros(uniq.max() + 1, np.int32)
        remap[uniq] = np.arange(1, uniq.size + 1)
        units = remap[units]
        n_units = int(uniq.size)
        del small_mask, iy_s, ix_s, uniq, remap
        log(f'    碎单元(<{UNIT_MIN_AREA}像元)已并入邻近, 保留 {n_units} 个')
    del small_ids
    gc.collect()

    u_sizes = np.bincount(units.ravel())
    u_area = u_sizes[1:] * (RES_C * RES_C) / 1e6
    log(f'    单元面积: 中位数 {np.median(u_area):.3f}km², p90 {np.percentile(u_area, 90):.2f}km², '
        f'最大 {u_area.max():.2f}km², >=5km²覆盖 {100*u_area[u_area>=5].sum()/u_area.sum():.1f}%')
    del u_sizes, u_area
    gc.collect()

    # 单元内高程归一化位置(平滑DEM上取min/max稳健化)
    ids = np.arange(1, n_units + 1)
    u_min = ndimage.minimum(dem_s, labels=units, index=ids)
    u_max = ndimage.maximum(dem_s, labels=units, index=ids)
    lut_min = np.zeros(n_units + 1, np.float64)
    lut_max = np.zeros(n_units + 1, np.float64)
    lut_min[1:] = u_min
    lut_max[1:] = u_max
    del u_min, u_max
    gc.collect()
    u_rng = lut_max[units] - lut_min[units]
    t25 = np.clip((dem_s - lut_min[units]) / np.maximum(u_rng, 1.0), 0.0, 1.0).astype(np.float32)
    del units, u_rng, lut_min, lut_max, dem_s
    gc.collect()
    log('    t_unit分位数: ' + str({p: round(float(np.percentile(t25, p)), 3)
                                   for p in [10, 25, 50, 75, 90]}))

    # 坝子TPI(物理值)
    tpi_basin25 = dem25 - ndimage.uniform_filter(dem25, int(BASIN_TPI_WIN / RES_C), mode='nearest')

    # ---------- 阶段2: 上采样至5m ----------
    log('[2] 上采样至5m...')
    t5 = upsample(t25, transform25)
    hand5 = upsample(hand25, transform25)
    del t25, hand25
    gc.collect()
    tpi_basin5 = upsample(tpi_basin25, transform25)
    relief2k5 = upsample(relief2k25, transform25)
    del tpi_basin25, relief2k25, dem25
    gc.collect()

    # ---------- 阶段3: 5m坡度 ----------
    log('[3] 5m坡度...')
    gy, gx = np.gradient(dem5, RES_F)
    slope5 = np.degrees(np.arctan(np.hypot(gx, gy))).astype(np.float32)
    del gy, gx
    gc.collect()

    # ---------- 阶段4: 坝区(开阔度第一性判据) ----------
    log('[4] 坝区识别...')
    basin = (slope5 < BASIN_SLOPE_TH) & (tpi_basin5 < BASIN_TPI_TH)
    del tpi_basin5
    gc.collect()
    basin = ndimage.binary_closing(basin, structure=np.ones((5, 5), bool))
    st8 = np.ones((3, 3), bool)
    # 面积过滤
    lab, n = ndimage.label(basin, structure=st8)
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep_ids = np.where(sizes >= BASIN_MIN_AREA / (RES_F * RES_F))[0]
    basin = np.isin(lab, keep_ids)
    del lab, sizes, keep_ids
    gc.collect()
    # 核心化: 逐像元半宽 >= BASIN_CORE_HALF_W (等价于半径该值的形态学开运算),
    # 剔除串珠细颈与窄谷带——坝子只留局部宽 >= 2*半宽 的开阔盆底
    dist_in = ndimage.distance_transform_edt(basin)
    core = basin & (dist_in >= BASIN_CORE_HALF_W / RES_F)
    del dist_in
    gc.collect()
    core = ndimage.binary_closing(core, structure=np.ones((5, 5), bool))
    lab, n = ndimage.label(core, structure=st8)
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep_ids = np.where(sizes >= BASIN_CORE_MIN_AREA / (RES_F * RES_F))[0]
    basin = np.isin(lab, keep_ids)
    n_core = int(keep_ids.size)
    del lab, sizes, keep_ids, core
    gc.collect()
    log(f'    坝区 {100*basin.mean():.2f}% (核心化[{BASIN_CORE_HALF_W:.0f}m半宽/宽>={2*BASIN_CORE_HALF_W:.0f}m]后 '
        f'{n_core} 个盆底)')

    # ---------- 阶段5: 精细亚类(5m产品) ----------
    log('[5] 地貌亚类栅格(5m)...')
    sub5 = np.zeros(dem5.shape, np.uint8)
    z_low = dem5 < HILL_Z_MAX
    sub5[z_low & (relief2k5 < RELIEF_LOW_HILL)] = 1
    sub5[z_low & (relief2k5 >= RELIEF_LOW_HILL)] = 2
    sub5[(dem5 >= HILL_Z_MAX) & (dem5 < 1000)] = 3
    sub5[(dem5 >= 1000) & (dem5 < 3500)] = 4
    sub5[(dem5 >= 3500) & (dem5 < 5000)] = 5
    sub5[dem5 >= 5000] = 6
    sub5[basin] = 7
    del relief2k5, z_low
    gc.collect()
    out_sub = os.path.join(OUT, 'geomorph_subclass_5m.tif')
    prof = dict(driver='GTiff', height=H5, width=W5, count=1, dtype='uint8',
                crs=crs5, transform=transform5, nodata=0, compress='lzw',
                tiled=True, blockxsize=512, blockysize=512)
    with rasterio.open(out_sub, 'w', **prof) as dst:
        dst.write(sub5, 1)
        dst.write_colormap(1, {c: SUB_COLORS[c] + (255,) for c in SUB_COLORS} | {0: (0, 0, 0, 0)})
    log(f'    {out_sub}')

    # ---------- 阶段6: 坡面位置参数 t 坡位判定 ----------
    log('[6] 坡位判定 (t=HAND/(HAND+RHAND))...')
    dn = t5 < POS_DN             # 坡下: 靠近排水线
    up = t5 > POS_UP             # 坡上: 靠近分水线
    dn[(slope5 < TERRACE_SLOPE) & (hand5 < TERRACE_HAND)] = True   # 阶地/坡麓归下部
    del t5, hand5, slope5
    gc.collect()
    mid = ~dn & ~up

    # ---------- 阶段7: 组合8类 ----------
    log('[7] 组合8类...')
    result = np.zeros(dem5.shape, np.uint8)
    result[basin] = 1
    hill_zone = dem5 < HILL_Z_MAX
    mtn_zone = ~hill_zone
    result[~basin & hill_zone & up] = 3
    result[~basin & hill_zone & mid] = 4
    result[~basin & hill_zone & dn] = 5
    result[~basin & mtn_zone & up] = 6
    result[~basin & mtn_zone & mid] = 7
    result[~basin & mtn_zone & dn] = 8
    del dn, up, mid, basin, hill_zone, mtn_zone, dem5
    gc.collect()

    # ---------- 阶段8: 后处理(坝子类保护: 不参与投票/去斑, 形态只由核心化决定) ----------
    log('[8] 众数滤波与小图斑处理(坝子类保护)...')
    slope_classes = [3, 4, 5, 6, 7, 8]
    # 坝子(1)不在 classes 中 -> in_domain=False -> 滤波自动保持坝子原值,
    # 且不参与邻近像元的多数投票, 消除滤波沿窄坡下条蔓延出伪坝子的通道
    result = majority_filter_3x3(result, slope_classes, 1)
    mask_basin = result == 1
    small = np.zeros(result.shape, bool)
    for c in slope_classes:
        m = result == c
        if not m.any():
            continue
        lab, _ = ndimage.label(m, structure=st8)
        sizes = np.bincount(lab.ravel())
        bad = np.where(sizes < 400)[0]
        bad = bad[bad != 0]
        if bad.size:
            small |= np.isin(lab, bad)
        del lab, sizes, bad, m
        gc.collect()
    if small.any():
        # 填充时把坝子区一并跳过: 小斑只从坡位类取值, 不会填充成坝子
        skip = small | mask_basin
        _, (iy, ix) = ndimage.distance_transform_edt(skip, return_indices=True)
        result[small] = result[iy[small], ix[small]]
        del iy, ix
    del small, mask_basin
    gc.collect()

    # ---------- 阶段9: 输出 ----------
    log('[9] 输出栅格...')
    out_tif = os.path.join(OUT, 'terrain_position_5m.tif')
    with rasterio.open(out_tif, 'w', **prof) as dst:
        dst.write(result, 1)
        dst.write_colormap(1, {c: CLASS_COLORS[c] + (255,) for c in CLASS_COLORS} | {0: (0, 0, 0, 0)})
    log(f'    {out_tif}')

    vals, cnts = np.unique(result[result > 0], return_counts=True)
    total = cnts.sum() * RES_F * RES_F
    lines = ['地形部位划分面积统计(编码沿用calc_dxbw.py, HAND/RHAND体系)', '=' * 52]
    for v, c in zip(vals, cnts):
        a = c * RES_F * RES_F
        lines.append(f'{v}  {CLASS_NAMES[v]:　<6} {a/1e6:>10.2f} km²  {100*a/total:6.2f}%')
    lines.append(f'合计 {total/1e6:.2f} km²')
    report = '\n'.join(lines)
    with open(os.path.join(OUT, 'class_report.txt'), 'w', encoding='utf-8') as f:
        f.write(report + '\n')
    log(report)

    np.save(os.path.join(WORK, 'result_5m.npy'), result)
    log('[完成] 栅格主流程结束。矢量化请运行 raster_to_shp.py')


if __name__ == '__main__':
    main()
