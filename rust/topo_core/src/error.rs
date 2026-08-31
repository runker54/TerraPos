use thiserror::Error;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("IO 错误: {0}")]
    Io(#[from] std::io::Error),
    #[error("GeoTIFF 错误: {0}")]
    Tiff(#[from] tiff::TiffError),
    #[error("{0}")]
    Invalid(String),
    #[error("处理已取消")]
    Cancelled,
}

pub type Result<T> = std::result::Result<T, CoreError>;
