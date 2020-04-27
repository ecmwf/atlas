Changelog
=========

All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.20.2] - 2019-04-27
### Fixed
- Avoid 100ds of compilation warnings introduced in version 0.20.0

## [0.20.1] - 2019-04-08
### Fixed
- Make feature BOUNDSCHECKING work again. It was not turned on for DEBUG builds
- Workaround clang OpenMP bug
- Fix Segfault due to unexpected order of destruction of singleton objects

### Added
- atlas-grids tool can now be used to compute approximate North-South grid resolution

## [0.20.0] - 2019-03-06
### Fixed
- Pole edges hould not be created for global regular grids with points at poles
- Update compatibility with more recent GridTools 
- HaloExchange with CUDA
- Self registration from external library

### Added
- Proj-based projections
- OpenMP functions for sorting, filling, copying
- Parallel structured grid interpolation

### Changed
- Grid iterators can have random access
- Speed improvements for StructuredColumns constructor
- Speed improvements for LonLatPolygon::contains()
- Port to ecbuild 3 (also minimum required version)
- Tidying up of atlas tools 

## [0.19.2] - 2020-01-28
### Changed
- Compatibility with eckit 1.7 due to API change in eckit::LocalConfiguration

## [0.19.1] - 2019-12-19
### Fixed
- Keep Gaussian identity of a Gaussian grid if a given domain does not crop any latitudes
- Fix naming for LegendreCache, to be more specific, and platform independent

## [0.19.0] - 2019-10-01
### Fixed
- Lambert ( conformal conic ) projection xy coordinates are now corrected

### Changed
- LambertProjection renamed to LambertConformalConic

### Added
- Reordering of nodes strategies (Hilbert curve, ReverseCuthillMckee)
- Preliminary CellColumns functionspace with Gmsh IO; halos are not yet fully supported


## [0.18.1] - 2019-08-10
### Fixed
- Match vertical structured interpolation to IFS
- Fix in creating vertical dimension in StructuredColumns using interval
- Fix in caching StructuredColumnsHaloExchange


## [0.18.0] - 2019-07-15
### Changed
- Make grid hashes crossplatform

### Added
- Fortran: Can now use TransLocal
- Fortran: Can now create unstructured grids
- LonLat Bounding box computations for grids using its projection
- Serialisation of Mesh Connectivity tables to/from eckit::Stream

### Fixed
- Structured interpolation bugs
- StructuredColumns bug with iSend
- Memory corruption in Spectral functionspace with GT CUDA backend


## [0.17.2] - 2019-06-04
### Fixed
- Compilation with PGI 19.4
- Structured Grid creation for periodic domains that do not start at 0 degrees longitude (Greenwich)


## [0.17.1] - 2019-04-22
### Added
- Option to declaration of field type (vector/scalar) when creating field in FunctionSpace
- New projection: Lambert Azimuthal Equal Area
- StructuredInterpolation2D to target FunctionSpace StructuredColumns

### Fixed
- Compilation with IBM XL 19
- Compilation with Intel 19


## [0.17.0] - 2019-04-02
### Changed
- OpenMP is now private dependency
- Dependencies are now added in a modern CMake3 way
- Fortran modules are installed in <install-prefix>/module/atlas

### Added
- Spectral functionspace better used in distributed context
- Nabla now holds shared_ptr to Method
- Fortran access to vertical coordinates from StructuredColumns

### Fixed
- Compilation with PGI/19.1 (regression)
- Add missing halo_exchange for Fortran rank-4 arrays
- Memory leaks with StructuredColumns

## [0.16.0] - 2019-02-14
### Changed
- Interpolation makes use of OpenMP
- Cleanup of header includes
- fypp Fortran preprocessor is ported to fckit 0.6

### Added
- Parallel structured interpolation methods (2D,3D): linear, cubic, quasicubic
- Interpolation for multi-level and multi-variable fields
- atlas_Trace: Fortran API and use within OpenMP parallel regions
- StructuredColumns halo-exchange for vector fields
- Field::halo_exchange() function

### Fixed
- Fortran compilation with PGI 18.10
- Access to Field view within OpenMP parallel region
- FunctionSpaces use only required halo, even if larger halo is available in mesh
- Fixed faulty name of a Field when created through Fortran API, wrapping existing memory
- Fix NodeColumns functionspace when mesh is created from projected grid.
- Parallel interpolation from regular lonlat grid.
- Spectral spherical harmonics transforms for large cases

## [0.15.2] - 2018-08-31
### Changed
- Initialisation of Fields to signalling NaN in debug builds, uninitialised in
  non-debug builds (used to be initialised to zero as part of std::vector construction)

### Added
- Implementation of cropped unstructured grids so that spectral transforms to
  unstructured grids are allowed

### Fixed
- Spectral transforms to grids including pole and equator
- Build with gridtools CUDA backend

## [0.15.1] - 2018-07-17
### Fixed
- Compilation for Intel 18 debug
- Memory bug for spherical harmonics
- Compatibility with fckit 0.5.1

## [0.15.0] - 2018-06-19
### Changed
- Native Array data storage uses now a raw C pointer instead of std::vector
- Significant performance improvements to Spherical harmonics transforms

### Fixed
- Various bugs related to parallel halos
- Bit reproducibility for parallel interpolation

## [0.14.0] - 2018-03-22
### Added
- Spherical Harmonics transforms can receive a cache memory handle

### Changed
- Earth interface (C++)
- Requires eckit 0.20.0, fckit 0.5.0

## [0.13.2] - 2018-03-20
### Fixed
- C++ compilation using PGI 17, 18 and GCC 7
- Support Python 3 to generate Fortran bindings
- Travis CI linked with github repository
- Problem with CUDA allocated memory

## [0.13.1] - 2018-03-01
### Fixed
- Fortran compilation using Intel 18
- GridTools compatibility

## 0.13.0 - 2018-02-16

[Unreleased]: https://github.com/ecmwf/atlas/compare/master...develop
[0.20.2]: https://github.com/ecmwf/atlas/compare/0.20.1...0.20.2
[0.20.1]: https://github.com/ecmwf/atlas/compare/0.20.0...0.20.1
[0.20.0]: https://github.com/ecmwf/atlas/compare/0.20.0...0.19.2
[0.19.2]: https://github.com/ecmwf/atlas/compare/0.19.1...0.19.2
[0.19.1]: https://github.com/ecmwf/atlas/compare/0.19.0...0.19.1
[0.19.0]: https://github.com/ecmwf/atlas/compare/0.18.1...0.19.0
[0.18.1]: https://github.com/ecmwf/atlas/compare/0.18.0...0.18.1
[0.18.0]: https://github.com/ecmwf/atlas/compare/0.17.2...0.18.0
[0.17.2]: https://github.com/ecmwf/atlas/compare/0.17.1...0.17.2
[0.17.1]: https://github.com/ecmwf/atlas/compare/0.17.0...0.17.1
[0.17.0]: https://github.com/ecmwf/atlas/compare/0.16.0...0.17.0
[0.16.0]: https://github.com/ecmwf/atlas/compare/0.15.2...0.16.0
[0.15.2]: https://github.com/ecmwf/atlas/compare/0.15.1...0.15.2
[0.15.1]: https://github.com/ecmwf/atlas/compare/0.15.0...0.15.1
[0.15.0]: https://github.com/ecmwf/atlas/compare/0.14.0...0.15.0
[0.14.0]: https://github.com/ecmwf/atlas/compare/0.13.2...0.14.0
[0.13.2]: https://github.com/ecmwf/atlas/compare/0.13.1...0.13.2
[0.13.1]: https://github.com/ecmwf/atlas/compare/0.13.0...0.13.1
