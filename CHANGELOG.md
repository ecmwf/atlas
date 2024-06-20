Changelog
=========

All notable changes to this project will be documented in this file.

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.38.0] - 2024-06-20

### Added
- Make non_linear interpolation independent of a chosen value type by @wdeconinck in https://github.com/ecmwf/atlas/pull/176
- Procedure to carry out a regridding from high to low resolution (binning) by @mo-lormi in https://github.com/ecmwf/atlas/pull/191
- Add Fortran interface for node-to-edge connectivity building by @benjaminmenetrier in https://github.com/ecmwf/atlas/pull/209
- CUDA/OpenACC capable fields with Native storage backend  @sbrdar and @wdeconinck in https://github.com/ecmwf/atlas/pull/182

### Changed
- Make non_linear interpolation independent of a chosen value type by @wdeconinck in https://github.com/ecmwf/atlas/pull/176

### Fixed
- Remove float in Triag2D intersection algorithm by @fmahebert in https://github.com/ecmwf/atlas/pull/203
- Allow zero-sized interpolation target functionspace by @odlomax in https://github.com/ecmwf/atlas/pull/206
- Made sure cubed-sphere interpolation method always sets metadata. by @odlomax in https://github.com/ecmwf/atlas/pull/208
- Avoid silent errors accessing Fieldset fields by ambiguous names by @wdeconinck in https://github.com/ecmwf/atlas/pull/210
- Fixes opposite pole coordinates by @benjaminmenetrier in https://github.com/ecmwf/atlas/pull/202

## [0.37.0] - 2024-04-09
### Added
- Add SphericalVector interpolation method using parallel transport (#163)
- Added arrayForEachDim method

### Fixed
- Bugfix for regional grids with ny > nx
- Fix build for configuration setting ATLAS_BITS_LOCAL=64

### Changed
- Use new LocalConfiguration baseclass functions in util::Config and util::Metadata instead of eckit::Value backdoor
- atlas_io is an adaptor library when eckit_codec is available (#181)

## [0.36.0] - 2023-12-11
### Added
- Add TriangularMeshBuilder with Fortran API, so far for serial meshes only
- Add Fortran API for CellCollumns functionspace
- Add EqualAreaPartitioner which is geometry based rather than loadbalanced like EqualRegionsPartitioner

### Fixed
- Compilation with Intel LLVM compiler
- Fix 180 degrees phase shift error in MDPI_gulfstream function

### Changed
- Update install scripts
- Preparation for using eckit::codec as backend for atlas_io

## [0.35.1] - 2023-10-24
### Added
- Don't output with output::Gmsh the triangle elements with wrong orientation when coordinates are "lonlat"
- Add control to skip Gmsh-output of triangles with too large edge length ratio
- Configurable geometry in KDTree
- Use configurable KDTree geometry in PointCloud
- atlas-grid-points executable to list grid points
- Allow constructor of atlas::io::ArrayShape with different integer types
- Support atlas_io with vector<std::int64_t>

### Fixed
- Control FPE in StructuredColumns::checksum to avoid overflow and invalid in unimportant halo regions
- Fixes to MeshBuilder validate_mesh_vs_grid

### Changed
- Use search radius in FiniteElement interpolation when mesh defines metadata to do so

## [0.35.0] - 2023-02-10
### Added
- Add accessors for the GridBox class (MIR-632)
- Add FunctionSpace::partition() field
- Add configurable MPI communicator to data structures (#131)
- Add single precision support for fvm::Nabla
- Implement PointCloud parallel construction with halo_radius

### Fixed
- Fix StructuredMeshGenerator for Gaussian grids that don't start at 0 degrees (ATLAS-375)
- Remove functionspace::Spectral via Trans initialisation (fixes #153)
- Enable MeshBuilder to set up the Grid (#152)
- Fix use of cmake option ATLAS_ENABLE_TRANS as in a bundle

### Changed
- Cleanup long-standing temporary element types. API change but with deprecated old API
- Deprecated rename Mesh::nb_partitions -> Mesh::nb_parts
- Implement Delaunay triangulation using Qhull instead of CGAL but CGAL can still be enabled instead for now

## [0.34.0] - 2023-07-10
### Added
- Fieldset::metadata (#126)
- Fieldset::adjointHaloExchange
- Field/Fieldset::clone method
- Functions to enable/disable FPE
- Add function to build mesh from imported connectivity data (#135)
- Implement field::for_each capabilities (#139)
- Introduce colon-separated environment variable ATLAS_PLUGIN_PATH to simplify plugin detection
- Introduce atlas::mdspan, contributed from github.com/kokkos/mdspan
- Add function Field::horizontal_dimension() -> std::vector<idx_t>
- Setup horizontal_dimensions() for BlockStructuredColumns fields
- Upgrade the halo exchange procedure for the function space 'PointCloud' (#120)

### Fixed
- Enable latitude normalisation in KDTree coordinate transform (#140)
- Fix LocalView indexing bug for non-contiguous slices
- C++17 flag public propagation to downstream C++ CMake packages
- Fix cells().global_index() metadata for RegularLonLat grids in StructuredMeshGenerator
- BuildHalo: mark interior added cells as ghost

## [0.33.0] - 2023-04-03
### Added
- Add support for StructuredPartitionPolygon with halo > 0
- Add information on atlas having PROJ support

### Changed
- C++17 standard is now a requirement

### Fixed
Fix StructuredInterpolation2D with retry for failed stencils

## [0.32.1] - 2023-02-09
### Added
- Added (lon, lat) to (alpha, beta) transforms to cubed sphere projection

## [0.32.0] - 2023-01-23
### Added
- Added BlockStructuredColumns FunctionSpace
- Added function SolidBodyRotation
- Added convenience Fortran constructors for ShiftedLonLat, SHiftedLon, ShiftedLat
- Support for more FunctionSpaces in various interpolation methods

### Changed
- PolygonLocator now wraps around longitudes for non-matching domains
- StructuredMeshGenerator can generate meshes with partitions without elements
- SerialDistribution makes every MPI task have the entire mesh

### Fixed
- Remove assertion for checking of empty mesh in NodeColumns FunctionSpace
- Gaussian latitudes can now be computed for very high resolution, without running out of memory.
- Interpolation to functionspaces with empty partitions

## [0.31.1] - 2022-11-11
### Fixed
- Fix bug introduced in StructuredMeshGenerator global numbering with release 0.31.0 in commit a63fc62a2
- Fix healpix global numbering for pentagon pole elements in parallel
- Fix validity check of atlas::HealpixGrid

## [0.31.0] - 2022-11-10
### Added
- Extend PointCloud functionspace to do halo exchanges, including Fortran API
- Add FunctionSpace::gather and FunctionSpace::scatter abstraction

### Changed
- Improve performance of MatchingMeshPartitionerLonLatPolygon using OpenMP
- Improve performance of BuildHalo using OpenMP and unordered_map
- Improve performance of StructuredMeshGenerator using OpenMP
- Improve performance of ATLAS_TRACE
- Reduce memory peak in GatherScatter setup
- Global element numbering of RegularLonLat grids is now following rows independent of partitioning

### Fixed
- Running CI with Github Actions
- Fix output of atlas-grids y-range for shifted grids
- Fix building with ATLAS_BITS_LOCAL=64


## [0.30.0] - 2022-08-22
### Added
- Fortran API for Interpolation::execute_adjoint()
- Fortran API for FieldSet::name()
- Fortran API for MeshGenerator::generate(grid,partitioner)
- Pentagon element type
- SphericalHarmonic function
- New interpolation method: ConservativeSphericalPolygonInterpolation
- Support 'variables' option in functionspace::PointCloud::createField

### Changed
- Atlas-IO is now standalone project, still embedded but only depending on eckit
- Deprecate Trans naming of 'ifs' or 'trans' in favour of 'ectrans'
- Default StructuredMeshGenerator partitioner is equal_regions instead of trans/ectrans

### Fixed
- Fix global numbering HEALPix grid to standard
- Fix NodeColumns remote_index for parallel orca grids
- Fix use of ATLAS_LINALG_DENSE_BACKEND environment variable

## [0.29.0] - 2022-04-21
### Added
- MatchingMeshPartitioner "cubedsphere"
- Interpolator "cubedsphere-bilinear"
- Improvements to Interpolation::Cache
- Add support for rank 2 fields when a nonlinear action is added to the interpolator
- Create Array using ArraySpec only

### Changed
- FieldSet::has(...) replaces FieldSet::has_field(...)
- Metadata return value to Interpolation::execute()
- Rename BilinearRemapping to UnstructuredBilinearLonLat

### Fixed
- Compatibility with proj version >= 8
- Compatibility with eckit version <= 1.18.5
- Compatibility with GridTools backend and using 64bit idx_t
- Wrongly computed Jacobian::transpose() introduced in 0.28.0
- Fix bug where using ectrans was not enabling adjoint of invtrans
- Avoid segfault when OpenMP tasking is broken, as it is with AppleClang and LLVM libomp

## [0.28.1] - 2022-03-14
### Fixed
- Fix compilation for GNU 7.3

## [0.28.0] - 2022-03-02
### Added
- Assignment of ArrayView from ArrayView
- Grid "regional_variable_resolution" via a new VariableResolutionProjection
- CubedSphereDualMeshGenerator
- VortexRollup function as analytical field for initialising data
- ConvexSphericalPolygon utility class
- Improve Projection::Jacobian
- Initial implementation for bilinear interpolation for unstructured meshes

### Changed
- Use new eckit (1.19.0) Sparse and Dense linear algebra API
- General robustness improvements to CubedSphere to using functionspaces with various halos

### Fixed
- Workarounds to fix compilation with Fujitsu compiler
- Workarounds to avoid Cray compiler problems with certain flag combinations
- CellColumns::haloExchange for meshes with multiple element types
- Computation of HEALPix mesh remote indices.


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
[0.38.0]: https://github.com/ecmwf/atlas/compare/0.37.0...0.38.0
[0.37.0]: https://github.com/ecmwf/atlas/compare/0.36.0...0.37.0
[0.36.0]: https://github.com/ecmwf/atlas/compare/0.35.1...0.36.0
[0.35.1]: https://github.com/ecmwf/atlas/compare/0.35.0...0.35.1
[0.35.0]: https://github.com/ecmwf/atlas/compare/0.34.0...0.35.0
[0.34.0]: https://github.com/ecmwf/atlas/compare/0.33.0...0.34.0
[0.33.0]: https://github.com/ecmwf/atlas/compare/0.32.1...0.33.0
[0.32.1]: https://github.com/ecmwf/atlas/compare/0.32.0...0.32.1
[0.32.0]: https://github.com/ecmwf/atlas/compare/0.31.1...0.32.0
[0.31.1]: https://github.com/ecmwf/atlas/compare/0.31.0...0.31.1
[0.31.0]: https://github.com/ecmwf/atlas/compare/0.30.0...0.31.0
[0.30.0]: https://github.com/ecmwf/atlas/compare/0.29.0...0.30.0
[0.29.0]: https://github.com/ecmwf/atlas/compare/0.28.1...0.29.0
[0.28.1]: https://github.com/ecmwf/atlas/compare/0.28.0...0.28.1
[0.28.0]: https://github.com/ecmwf/atlas/compare/0.27.0...0.28.0
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
