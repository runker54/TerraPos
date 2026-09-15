# -*- coding: utf-8 -*-
"""合成 DEM 场景生成器(验收体系的数据源)。

所有常量在文件顶部; 不使用命令行参数。
运行: D:/worker_code/.venvgis/Scripts/python.exe python/build_synthetic_dem.py

场景均为解析定义的确定性几何(无随机成分; 噪声变体使用固定种子),
每个场景输出:
  <名>.tif            float32 DEM(投影米制, CGCS2000 3度带 EPSG:4545 仿射)
  <名>_expected.tif   byte 期望语义区(1 盆地 2 上部 3 中部 4 下部 0 未定义)
"""
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

OUTPUT_DIR = Path(__file__).resolve().parents[1] / "data" / "synthetic" / "generated"
NODATA = -9999.0
EPSG = 4545  # CGCS2000 / 3-degree Gauss-Kruger CM 108E (投影, 米)
NOISE_SEED = 42
NOISE_SIGMAS = (0.25, 0.5)
RESOLUTIONS_M = (10.0, 25.0)  # 尺度稳定性采样的分辨率组


def _writer(path: Path, arr, res, expected=False):
    """写出 GeoTIFF(投影米制仿射, 原点任意但固定)。"""
    rows, cols = arr.shape
    driver = gdal.GetDriverByName("GTiff")
    dtype = gdal.GDT_Byte if expected else gdal.GDT_Float32
    ds = driver.Create(str(path), cols, rows, 1, dtype)
    ds.SetGeoTransform((500000.0, res, 0.0, 3000000.0, 0.0, -res))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(EPSG)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    if expected:
        band.SetNoDataValue(0)
    else:
        band.SetNoDataValue(NODATA)
    band.WriteArray(arr)
    band.FlushCache()
    ds = None


def _grid(size, res):
    y, x = np.mgrid[0:size, 0:size]
    return (x + 0.5) * res, (y + 0.5) * res


def broad_basin(res=10.0, size=401, shift=0.0):
    x, y = _grid(size, res)
    cx = size * res * 0.5
    r = np.hypot(x - cx, y - cx)
    floor = 800.0 + 0.0005 * x
    rim = np.clip((r - 700.0) / 250.0, 0.0, 1.0)
    z = floor + 120.0 * rim * rim + shift
    expected = np.zeros_like(z, dtype=np.uint8)
    expected[r < 400.0] = 1
    expected[(r >= 400.0) & (r < 650.0)] = 4
    expected[r >= 900.0] = 2
    expected[(r >= 650.0) & (r < 900.0)] = 3
    return z, expected


def v_valley(res=10.0, size=201, mirror=False):
    x, y = _grid(size, res)
    cx = size * res * 0.5
    xi = x if not mirror else (size * res - x)
    d = np.abs(xi - cx)
    z = 800.0 + 0.10 * d + 0.004 * y
    expected = np.zeros_like(z, dtype=np.uint8)
    expected[d < size * res * 0.1] = 4
    expected[d >= size * res * 0.35] = 2
    expected[(d >= size * res * 0.1) & (d < size * res * 0.35)] = 3
    return z, expected


def nested_ridges(res=10.0, size=257):
    x, y = _grid(size, res)
    ridge = 40.0 * (1.0 - np.cos(2.0 * np.pi * x / 640.0))
    z = 850.0 + ridge + 0.008 * y
    ridge_mask = ridge > 30.0
    valley_mask = ridge < 8.0
    expected = np.zeros_like(z, dtype=np.uint8)
    expected[valley_mask] = 4
    expected[ridge_mask] = 2
    expected[~(ridge_mask | valley_mask)] = 3
    return z, expected


def conical_hill(res=10.0, size=301):
    x, y = _grid(size, res)
    cx = size * res * 0.5
    r = np.hypot(x - cx, y - cx)
    z = 700.0 + 150.0 * np.maximum(0.0, 1.0 - r / 600.0)
    expected = np.zeros_like(z, dtype=np.uint8)
    expected[r < 120.0] = 2
    expected[(r >= 120.0) & (r < 480.0)] = 3
    expected[r >= 480.0] = 4
    return z, expected


def closed_pit(res=10.0, size=201):
    x, y = _grid(size, res)
    cx = size * res * 0.5
    r = np.hypot(x - cx, y - cx)
    z = np.where(
        r < 400.0,
        810.0 + 2.0 * (r / 400.0) ** 2,
        np.where(r < 800.0, 812.0 + 60.0 * ((r - 400.0) / 400.0) ** 2, 872.0 + 0.08 * (r - 800.0)),
    )
    expected = np.zeros_like(z, dtype=np.uint8)
    expected[r < 350.0] = 1
    expected[r >= 850.0] = 2
    return z, expected


def build_all():
    """生成全部场景与变体。"""
    if OUTPUT_DIR.exists():
        for p in OUTPUT_DIR.iterdir():
            p.unlink()
    else:
        OUTPUT_DIR.mkdir(parents=True)
    assert str(OUTPUT_DIR).replace("\\", "/").endswith("data/synthetic/generated")

    scenes = {}

    z, e = broad_basin()
    scenes["broad_basin"] = (z, e, 10.0)
    scenes["broad_basin_shift50"] = (z + 50.0, e, 10.0)          # 垂直平移
    z25, e25 = broad_basin(res=25.0, size=161)
    scenes["broad_basin_25m"] = (z25, e25, 25.0)                  # 分辨率稳定性

    z, e = v_valley()
    scenes["v_valley"] = (z, e, 10.0)
    zm, em = v_valley(mirror=True)
    scenes["v_valley_mirror"] = (zm, em, 10.0)                    # 镜像

    z, e = nested_ridges()
    scenes["nested_ridges"] = (z, e, 10.0)
    scenes["nested_ridges_rot90"] = (np.rot90(z, 1).copy(), np.rot90(e, 1).copy(), 10.0)  # 旋转
    rng = np.random.default_rng(NOISE_SEED)
    for sigma in NOISE_SIGMAS:
        zn = z + rng.normal(0.0, sigma, z.shape)
        scenes[f"nested_ridges_noise{sigma}"] = (zn, e, 10.0)     # 噪声
    z25, e25 = nested_ridges(res=25.0, size=103)
    scenes["nested_ridges_25m"] = (z25, e25, 25.0)

    z, e = conical_hill()
    scenes["conical_hill"] = (z, e, 10.0)
    z, e = closed_pit()
    scenes["closed_pit"] = (z, e, 10.0)

    zm, em = v_valley(res=10.0, size=129)
    scenes["ref_small"] = (zm, em, 10.0)  # 参考计算小样区

    for name, (z, e, res) in scenes.items():
        _writer(OUTPUT_DIR / f"{name}.tif", z.astype(np.float32), res)
        _writer(OUTPUT_DIR / f"{name}_expected.tif", e, res, expected=True)
    print(f"生成 {len(scenes)} 个场景 -> {OUTPUT_DIR}")


if __name__ == "__main__":
    gdal.UseExceptions()
    build_all()
