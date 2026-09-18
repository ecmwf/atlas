# atlas-sys

Build-only Rust crate for ECMWF's [atlas](https://github.com/ecmwf/atlas) (grid/mesh) C++ library.

This crate has **no Rust API**. It builds (or locates) the atlas C++ library and exports `DEP_ATLAS_SYS_ROOT` / `DEP_ATLAS_SYS_INCLUDE` for downstream `-sys` crates that need to link against atlas.

## Features

### Build strategy (mutually exclusive)

- `vendored` - Build atlas (and its eckit dependency) from source.
- `system` - Link against system-installed atlas.

`vendored` is enabled by default.

### Optional

All off by default, and only meaningful for `vendored` builds - a `system` build
gets whatever the installed atlas was compiled with. Enabling one whose library
CMake cannot find fails the build rather than silently dropping the feature.

- `omp` - OpenMP support. Requires `libomp` (e.g. `libomp-dev` on Ubuntu).
- `tesselation` - Unstructured mesh generation. Requires Qhull.
- `eigen` - Eigen linear algebra backend. Requires Eigen3.
- `fftw` - FFTW backend for spectral transforms. Requires FFTW with `double`.

## License

Apache-2.0
