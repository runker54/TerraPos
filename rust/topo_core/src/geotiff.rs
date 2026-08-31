//! 轻量 GeoTIFF 读写器（无 GDAL 依赖）
//! 读: 自研 IFD 解析(地理标签) + tiff crate 像素解码(LZW/Deflate/PackBits)
//! 写: 自研条带式 Deflate GeoTIFF (float32 / uint8+色表)

use crate::error::{CoreError, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom};
use std::path::Path;

/// 地理元信息（从源 DEM 提取并透传到输出）
#[derive(Debug, Clone, Default)]
pub struct GeoMeta {
    pub width: u32,
    pub height: u32,
    pub pixel_scale: [f64; 3],
    pub tiepoint: [f64; 6],
    pub geo_keys: Vec<u16>,
    pub geo_ascii: Option<String>,
}

impl GeoMeta {
    /// 由仿射参数构造（左上角 x/y + 分辨率）
    pub fn from_origin(width: u32, height: u32, x0: f64, y0: f64, res: f64) -> Self {
        GeoMeta {
            width,
            height,
            pixel_scale: [res, res, 0.0],
            tiepoint: [0.0, 0.0, 0.0, x0, y0, 0.0],
            geo_keys: Vec::new(),
            geo_ascii: None,
        }
    }
    pub fn origin(&self) -> (f64, f64) {
        (self.tiepoint[3], self.tiepoint[4])
    }
    pub fn resolution(&self) -> f64 {
        self.pixel_scale[0]
    }
}

#[derive(Debug, Clone)]
struct IfdEntry {
    typ: u16,
    #[allow(dead_code)]
    count: u64,
    data: Vec<u8>,
}

const TYPE_SIZE: [usize; 13] = [0, 1, 1, 2, 4, 8, 8, 1, 1, 2, 4, 8, 8];

/// 解析 TIFF IFD 全部条目（仅支持小端 II 格式，GeoTIFF 通用）
fn parse_ifd(f: &mut File) -> Result<HashMap<u16, IfdEntry>> {
    let mut head = [0u8; 8];
    f.read_exact(&mut head)?;
    if &head[0..2] != b"II" {
        return Err(CoreError::Invalid("仅支持小端(II) TIFF".into()));
    }
    let magic = u16::from_le_bytes([head[2], head[3]]);
    if magic != 42 {
        return Err(CoreError::Invalid("非 TIFF 文件".into()));
    }
    let ifd_off = u32::from_le_bytes(head[4..8].try_into().unwrap()) as u64;
    f.seek(SeekFrom::Start(ifd_off))?;
    let mut cnt_b = [0u8; 2];
    f.read_exact(&mut cnt_b)?;
    let cnt = u16::from_le_bytes(cnt_b) as u64;
    let mut map = HashMap::new();
    for _ in 0..cnt {
        let mut e = [0u8; 12];
        f.read_exact(&mut e)?;
        let tag = u16::from_le_bytes(e[0..2].try_into().unwrap());
        let typ = u16::from_le_bytes(e[2..4].try_into().unwrap());
        let count = u32::from_le_bytes(e[4..8].try_into().unwrap()) as u64;
        if !(1..13).contains(&typ) {
            continue;
        }
        let total = TYPE_SIZE[typ as usize] * count as usize;
        let data = if total <= 4 {
            e[8..8 + total].to_vec()
        } else {
            let off = u32::from_le_bytes(e[8..12].try_into().unwrap()) as u64;
            let pos = f.stream_position()?;
            f.seek(SeekFrom::Start(off))?;
            let mut buf = vec![0u8; total];
            f.read_exact(&mut buf)?;
            f.seek(SeekFrom::Start(pos))?;
            buf
        };
        map.insert(tag, IfdEntry { typ, count, data });
    }
    Ok(map)
}

fn entry_doubles(e: &IfdEntry) -> Vec<f64> {
    e.data
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn entry_u16s(e: &IfdEntry) -> Vec<u16> {
    e.data
        .chunks_exact(2)
        .map(|c| u16::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

/// 读取 float32 GeoTIFF（含地理元信息）
pub fn read_f32<P: AsRef<Path>>(path: P) -> Result<(Vec<f32>, GeoMeta)> {
    let path = path.as_ref();
    let mut f = File::open(path)?;
    let tags = parse_ifd(&mut f)?;
    let get_dim = |tag: u16, name: &str| -> Result<u32> {
        let e = tags
            .get(&tag)
            .ok_or_else(|| CoreError::Invalid(format!("缺 {}", name)))?;
        if e.typ == 3 {
            Ok(entry_u16s(e)[0] as u32)
        } else {
            Ok(e.data
                .chunks_exact(4)
                .map(|c| u32::from_le_bytes(c.try_into().unwrap()))
                .next()
                .unwrap_or(0))
        }
    };
    let width = get_dim(256, "ImageWidth")?;
    let height = get_dim(257, "ImageLength")?;
    let bits = tags
        .get(&258)
        .map(entry_u16s)
        .ok_or_else(|| CoreError::Invalid("缺 BitsPerSample".into()))?;
    let sfmt = tags.get(&339).map(entry_u16s).unwrap_or_default();
    if bits.first() != Some(&32) || sfmt.first() != Some(&3) {
        return Err(CoreError::Invalid(format!(
            "仅支持 float32 GeoTIFF (bits={:?}, sample_format={:?})",
            bits, sfmt
        )));
    }
    let mut meta = GeoMeta { width, height, ..Default::default() };
    if let Some(e) = tags.get(&33550) {
        let d = entry_doubles(e);
        if d.len() >= 3 {
            meta.pixel_scale = [d[0], d[1], d[2]];
        }
    }
    if let Some(e) = tags.get(&33922) {
        let d = entry_doubles(e);
        if d.len() >= 6 {
            meta.tiepoint = [d[0], d[1], d[2], d[3], d[4], d[5]];
        }
    }
    if let Some(e) = tags.get(&34735) {
        meta.geo_keys = entry_u16s(e);
    }
    if let Some(e) = tags.get(&34737) {
        meta.geo_ascii =
            Some(String::from_utf8_lossy(&e.data).trim_end_matches('\0').to_string());
    }

    let file = File::open(path)?;
    let mut dec = tiff::decoder::Decoder::new(BufReader::new(file))?
        .with_limits(tiff::decoder::Limits::unlimited());
    match dec.read_image()? {
        tiff::decoder::DecodingResult::F32(buf) => Ok((buf, meta)),
        other => Err(CoreError::Invalid(format!("非 float32 像素: {:?}", other))),
    }
}

// ---------------- 写出 ----------------
// 像素编码交由 tiff crate(zlib 压缩), GeoTags 通过 DirectoryEncoder::write_tag 注入

use tiff::encoder::colortype;
use tiff::tags::Tag;

/// 写 GeoTIFF: float32 高程栅格
pub fn write_f32<P: AsRef<Path>>(
    path: P,
    meta: &GeoMeta,
    data: &[f32],
) -> Result<()> {
    let file = File::create(path)?;
    let mut enc = tiff::encoder::TiffEncoder::new(BufWriter::new(file))?;
    {
        let mut img = enc.new_image::<colortype::Gray32Float>(meta.width, meta.height)?;
        let dir = img.encoder();
        dir.write_tag(Tag::Unknown(33550), &meta.pixel_scale[..])?;
        dir.write_tag(Tag::Unknown(33922), &meta.tiepoint[..])?;
        if !meta.geo_keys.is_empty() {
            dir.write_tag(Tag::Unknown(34735), &meta.geo_keys[..])?;
        }
        if let Some(a) = &meta.geo_ascii {
            dir.write_tag(Tag::Unknown(34737), a.as_str())?;
        }
        img.write_data(data)?;
    }
    Ok(())
}

/// 写 GeoTIFF: uint8 分类栅格(带 256 色表)
pub fn write_u8_cmap<P: AsRef<Path>>(
    path: P,
    meta: &GeoMeta,
    data: &[u8],
    cmap: &[[u8; 3]; 256],
) -> Result<()> {
    let file = File::create(path)?;
    let mut enc = tiff::encoder::TiffEncoder::new(BufWriter::new(file))?;
    {
        let mut img = enc.new_image::<colortype::Gray8>(meta.width, meta.height)?;
        let dir = img.encoder();
        // 色表: tag 320, SHORT[768], 通道优先 R,G,B, 值域 0..65535
        let mut cm: Vec<u16> = Vec::with_capacity(768);
        for ch in 0..3 {
            for px in cmap {
                cm.push((px[ch] as u16) * 257);
            }
        }
        dir.write_tag(Tag::Unknown(320), &cm[..])?;
        dir.write_tag(Tag::Unknown(33550), &meta.pixel_scale[..])?;
        dir.write_tag(Tag::Unknown(33922), &meta.tiepoint[..])?;
        if !meta.geo_keys.is_empty() {
            dir.write_tag(Tag::Unknown(34735), &meta.geo_keys[..])?;
        }
        if let Some(a) = &meta.geo_ascii {
            dir.write_tag(Tag::Unknown(34737), a.as_str())?;
        }
        img.write_data(data)?;
    }
    Ok(())
}


/// 预览辅助: 读取 uint8 分类栅格并盒式下采样到指定宽度
pub fn read_u8_preview<P: AsRef<Path>>(path: P, max_w: usize) -> Result<(Vec<u8>, usize, usize)> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let mut dec = tiff::decoder::Decoder::new(BufReader::new(file))?
        .with_limits(tiff::decoder::Limits::unlimited());
    let (w, h) = dec.dimensions()?;
    let data = match dec.read_image()? {
        tiff::decoder::DecodingResult::U8(buf) => buf,
        other => return Err(CoreError::Invalid(format!("非 uint8 像素: {:?}", other))),
    };
    let (w, h) = (w as usize, h as usize);
    if w <= max_w {
        return Ok((data, w, h));
    }
    let step = w / max_w + 1;
    let nw = w / step;
    let nh = h / step;
    let mut out = vec![0u8; nw * nh];
    for y in 0..nh {
        for x in 0..nw {
            out[y * nw + x] = data[(y * step) * w + x * step];
        }
    }
    Ok((out, nw, nh))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_f32() {
        let dir = std::env::temp_dir().join("topo_core_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("rt_f32.tif");
        let meta = GeoMeta::from_origin(64, 48, 371242.5, 3091607.5, 25.0);
        let data: Vec<f32> = (0..64 * 48).map(|i| (i % 97) as f32).collect();
        write_f32(&p, &meta, &data).unwrap();
        let (back, m2) = read_f32(&p).unwrap();
        assert_eq!(back.len(), data.len());
        assert_eq!(m2.width, 64);
        assert_eq!(m2.resolution(), 25.0);
        assert_eq!(m2.origin(), (371242.5, 3091607.5));
        assert!(back.iter().zip(data.iter()).all(|(a, b)| a == b));
    }
}
