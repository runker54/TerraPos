//! topo_core: 地形部位划分核心算法库（纯 Rust，无 GDAL 依赖）
//!
//! 模块:
//! - [`geotiff`]：轻量 GeoTIFF 读写
//! - [`filter`]：并行分离滤波(均值/min/max)
//! - [`terrain`]：坡度/起伏度/TPI
//! - [`hydro`]：Priority-Flood 填洼 + D8 流向/汇流
//! - [`distance`]：精确欧氏距离变换(含最近源索引)
//! - [`segment`]：峰顶提取 + 分水岭分割 + 连通域
//! - [`classify`]：地貌亚类/地形部位判据
//! - [`post`]：众数滤波/小图斑处理
//! - [`pipeline`]：参数模型 + 全流程编排

pub mod distance;
pub mod error;
pub mod filter;
pub mod geotiff;
pub mod hydro;
pub mod pipeline;
pub mod segment;
pub mod terrain;
