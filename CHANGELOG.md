Changelog
=========

All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.27.0] - 2021-12-03
### Added
- Adjoint interpolation with some restrictions
- Cubed sphere grid partitioner
- Cubed sphere parallel mesh generation
- Cubed sphere function spaces
- Fortran interfaces to Projection methods
- Support discovery of open-source ectrans
- Dense linear Algebra matrix_multiply abstraction

### Changed
- Remove etc/atlas/config.yaml because defaults should be in code
- Naming of sparse_matrix_multiply backend 'omp' -> 'openmp'
- Applied clang-format 13.0.0 (all files touched)

## [0.26.0] - 2021-08-23
### Added
- Support for cubed sphere grids and preliminary support for cubes sphere mesh generation.

### Fixed
- Compilation with altas_bits_local=64
- Too aggressive optimisation with gnu 11
- Compatibility with cmake 3.20 and nvhpc compilers


## [0.25.0] - 2021-05-18
### Added
- New concept Redistribution to redistribute field on same grid but with different partitioner
- Support for Trans::invtrans_adj() and related methods
- Introduce ~atlas/etc/atlas/config.yaml
- atlas-io with support for encoding/decoding std::array
- atlas-io print support for tiny arrays

### Changed
- atlas-io version 0.2
- Move util::SpecRegistry<T> to non-templated grid::SpecRegistry

### Fixed
- minor bug fixes

## [0.24.1] - 2021-04-06
### Fixed
- Periodic halo for HEALPix mesh
- General warnings and suggestions by DeepCode bot
- Compilation problems with clang 3.8


## [0.24.0] - 2021-03-19
### Fixed
- Fixed hang in band distribution for L160000x8000 grid
- Fixes to Spectral functionspace and TransIFS regarding vertical levels

### Changed
- Requires eckit 1.16

### Added
- atlas-io first version, not for operational use
- Spec registry for grids

## [0.23.0] - 2021-01-19
### Fixed
- Structured interpolation method interpolating to area straddling Greenwich.
- Fixes when compiling with ATLAS_BITS_LOCAL=64

### Changed
- Possibility to link to alternative open-source version of IFS trans library.

### Added
- Caching mechanism for interpolation

## [0.22.1] - 2020-10-22
### Fixed
- Installation of PGI compilers via tools/install-pgi.sh
- Allow dependency on older Eigen 3.2 which does not use CMake targets

## [0.22.0] - 2020-10-14
### Fixed
- Feature INIT_SNAN was not working
- Support KNearestNeighbour interpolation for functionspace with
  smaller halo than the mesh contains
- Support array size up to size_t limit

### Changed
- Migration to use ecbuild 3.4
- ATLAS_BITS_LOCAL can be configured to 32 or 64

### Added
- Fields can be created with given alignment, which adds padding in innermost dimension
- Added conservative interpolation with "grid-box average" and "grid-box maximum"
- Missing value definition for fields
- Add support for missing values in matrix-based interpolation methods
- Floating point trapping and signal handling mechanism
- Fortran: GridDistribution constructors
- Fortran: Domain access
- Fortran: Get lonlat_bounding_box via domain
- Possibility to access Jacobian of projections (with only some projections implemented)

## [0.21.0] - 2020-06-23
### Fixed
- Fixed Rotation order of applying the rotation angle
- Fix nabla computation of surface vector fields
- Fix registration and destruction issues of halo-exchange caches
- Workaround Clang compiler problem with OpenMP tasking, using introspection
- Fix bug in conversion of negative degrees to microdegrees
- Fix problem in distribution of work amongst OpenMP threads in StructuredColumns::setup
- Fix problem with StructuredColumns creation for grids defined with domains with negative West
- Fix computation of Grid::lonlatBoundingBox for arbitrary projections crossing the dateline.

### Changed
- Snap LinearSpacing values to start and endpoint to allow exact comparisons
- Improved performance and memory requirement of cropping of large StructuredGrids
- Regional grids by default now have a positive y-numbering (previously negative).

### Added
- PolygonLocator based on functionspace partition polygons.
- KDTree class which abstracts eckit concrete implementations, including Fortran interface
- Geometry class with factory mechanism to configure which spheroid to use, including Fortran interface
- StrcutredGrid function for i,j to index and inverse
- Fortran interface to create StructuredColumns with custom distribution
- Fortran interface to grid spec
- Fortran interface to Projection
- Adjoint of HaloExchange
- Plugin mechanism to load plugins at runtime.
- Fortran constructors for atlas_RegionalGrid
- Fortran constructors for projected reduced Gaussian grids
- Add copy constructor and assignment operator for atlas::vector
- Mercator projection support for scaling, and operation on ellipsoid.
- Grid Distribution can now also be created as a function, e.g. for Serial or Bands



## [0.20.2] - 2020-04-27
### Fixed
- Avoid 100ds of compilation warnings introduced in version 0.20.0

## [0.20.1] - 2020-04-08
### Fixed
- Make feature BOUNDSCHECKING work again. It was not turned on for DEBUG builds
- Workaround clang OpenMP bug
- Fix Segfault due to unexpected order of destruction of singleton objects

### Added
- atlas-grids tool can now be used to compute approximate North-South grid resolution

## [0.20.0] - 2020-03-06
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
[0.27.0]: https://github.com/ecmwf/atlas/compare/0.26.0...0.27.0
[0.26.0]: https://github.com/ecmwf/atlas/compare/0.25.0...0.26.0
[0.25.0]: https://github.com/ecmwf/atlas/compare/0.24.1...0.25.0
[0.24.1]: https://github.com/ecmwf/atlas/compare/0.24.0...0.24.1
[0.24.0]: https://github.com/ecmwf/atlas/compare/0.23.0...0.24.0
[0.23.0]: https://github.com/ecmwf/atlas/compare/0.22.1...0.23.0
[0.22.1]: https://github.com/ecmwf/atlas/compare/0.22.0...0.22.1
[0.22.0]: https://github.com/ecmwf/atlas/compare/0.21.0...0.22.0
[0.21.0]: https://github.com/ecmwf/atlas/compare/0.20.2...0.21.0
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
