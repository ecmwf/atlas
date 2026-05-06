# atlas-sys

Build-only Rust crate for ECMWF's [atlas](https://github.com/ecmwf/atlas) (grid/mesh) C++ library.

This crate has **no Rust API**. It builds (or locates) the atlas C++ library and exports `DEP_ATLAS_SYS_ROOT` / `DEP_ATLAS_SYS_INCLUDE` for downstream `-sys` crates that need to link against atlas.

## Features

### Build strategy (mutually exclusive)

- `vendored` - Build atlas (and its eckit dependency) from source.
- `system` - Link against system-installed atlas.

`vendored` is enabled by default.

## License

Apache-2.0
