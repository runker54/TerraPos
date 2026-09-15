# TerraPos 地形部位划分工具

[简体中文](README.md) | [English](README_EN.md)

只依赖一份投影米制 DEM，自适应划分中国西南复杂山地的地形部位，
输出 7 类有效结果与全部诊断图层。

## 地形部位编码

| 编码 | 部位 | 语义 |
|---|---|---|
| 0 | NoData | 原 DEM 无效区，全程保留 |
| 1 | 山间/宽谷盆地 | 盆地对象（宽度/面积/内起伏/围限/水文连通多重检验） |
| 2 | 保留编码 | 永不生成 |
| 3/4/5 | 丘陵上部/中部/下部 | 海拔 <500m，模糊隶属度上/中/下 |
| 6/7/8 | 山地坡上/坡中/坡下 | 海拔 ≥500m，模糊隶属度上/中/下 |

地貌亚类栅格独立输出：低丘 / 高丘 / 低山 / 中山 / 高山 / 极高山 / 盆地。

## 使用

```bash
cd rust
cargo build --release
./target/release/topo_app.exe
```

界面：导入 DEM → 选择分析方案（可默认）→ 运行 → 预览结果与诊断图层。
输出目录包含：

- `terrain_position.tif` 最终 7 类部位
- `geomorph_subclass.tif` 地貌亚类
- `terrain_confidence.tif` 置信度(0-100)
- `class_report.txt` 面积/对象/置信度/阶段耗时报告
- `diagnostics/` 11 个诊断图层(可在界面关闭)

## 算法要点

1. 输入校验(投影米制) + 边界 NoData 保留 + <=0.01km2 内部孔填补；
2. 地貌/水文双表面：保形平滑 + 有界 Priority-Flood(z-limit, 深洼保留)；
3. 米制尺度族(125-4000m)稳健金字塔 -> 特征尺度自适应选择；
4. 四级嵌套河网(物理汇水面积 0.05/0.2/1/5 km2) + 子流域边界脊线；
5. 坡面单元内约束距离 dv/dr -> 相对位置 q + 低起伏权重；
6. geomorphon 十标准形态 x 规格证据表 -> 模糊上/中/下隶属度；
7. 盆地对象：低平候选 -> 宽度核心 -> 多重检验 -> 受约束完整边界重建；
8. 受约束清理：单元内小斑归并 + 脊->谷轨迹单调 DP 修正。

## 验证

`python/` 目录提供合成场景、参考公式与验收脚本(见 python/README.md)；
无人工真值时不声称分类精度，验收依赖解析预期、拓扑一致性与尺度稳定性。

## 发布

发布包（TerraPos-v0.0.1-win64.zip）见 GitHub [Releases](https://github.com/runker54/TerraPos/releases)，解压即用；
发布说明见 docs/RELEASE_NOTES-v0.0.1.md。
