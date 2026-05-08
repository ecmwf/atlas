Atlas Environment Variables {#atlas_environment_variables}
===========================

This file lists Atlas environment variables found in code under atlas/src (and test-only ones under atlas/src/tests).

Notes:
- Boolean values are parsed by eckit translators. In practice use one of: 0/1, false/true, off/on, no/yes.
- List values use eckit vector parsing (typically comma-separated), for example: 10,20,30.

Core Logging And Runtime
------------------------

ATLAS_INFO
: Type: bool
: Default: true
: Effect: Enables Atlas info logging channel.
: Values:
  - true/1/on/yes: info logs enabled.
  - false/0/off/no: info logs disabled.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_WARNING
: Type: bool
: Default: true
: Effect: Enables Atlas warning logging channel.
: Values:
  - true/1/on/yes: warning logs enabled.
  - false/0/off/no: warning logs disabled.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_TRACE
: Type: bool
: Default: false
: Effect: Enables Atlas trace logging channel.
: Values:
  - true/1/on/yes: trace logs enabled.
  - false/0/off/no: trace logs disabled.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_TRACE_MEMORY
: Type: bool
: Default: false
: Effect: Enables memory tracing/reporting path in Atlas tracing.
: Values:
  - true/1/on/yes: memory tracing enabled.
  - false/0/off/no: memory tracing disabled.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_TRACE_BARRIERS
: Type: bool
: Default: false
: Effect: Enables barrier tracing.
: Values:
  - true/1/on/yes: barrier tracing enabled.
  - false/0/off/no: barrier tracing disabled.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_TRACE_REPORT
: Type: bool
: Default: false
: Effect: Prints trace report during finalise() when tracing is compiled in.
: Values:
  - true/1/on/yes: emit final trace report.
  - false/0/off/no: do not emit report.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_LOG_RANK
: Type: int
: Default: 0
: Effect: MPI rank that keeps normal console logging in several code paths.
: Values:
  - Any integer rank.
: Source Location: atlas/src/atlas/library/Library.cc, atlas/src/atlas/library/FloatingPointExceptions.cc, atlas/src/atlas/runtime/AtlasTool.cc

ATLAS_LOG_FILE
: Type: bool
: Default: false
: Effect: In AtlasTool-based programs, write rank-specific log files.
: Values:
  - true/1/on/yes: write log files (name derived from executable/display name and rank).
  - false/0/off/no: keep logging policy without per-rank logfile output.
: Source Location: atlas/src/atlas/runtime/AtlasTool.cc

ATLAS_WORKDIR
: Type: path string
: Default: current working directory
: Effect: AtlasTool logfile directory base.
: Values:
  - Any filesystem path.
: Source Location: atlas/src/atlas/runtime/AtlasTool.cc

ATLAS_FINALISES_MPI
: Type: bool
: Default: false
: Effect: If enabled, atlas::Library::finalise() calls atlas::mpi::finalise().
: Values:
  - true/1/on/yes: Atlas finalise also finalises MPI.
  - false/0/off/no: Atlas does not finalise MPI here.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_PLUGIN_PATH
: Type: string
: Default: empty
: Effect: Extra plugin search path(s) for eckit plugin manager (if supported by linked eckit).
: Values:
  - Empty: no extra plugin path.
  - Path list (colon-separated on Unix/macOS).
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_DATA_PATH
: Type: path list string
: Default: empty (Atlas also adds built-in fallback path entries)
: Effect: Additional data search path entries.
: Values:
  - Path list (colon-separated on Unix/macOS).
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_CACHE_PATH
: Type: path string
: Default: /tmp/cache
: Effect: Cache location used by atlas::Library::cachePath().
: Values:
  - Any filesystem path.
: Source Location: atlas/src/atlas/library/Library.cc

Floating Point And Signals
--------------------------

ATLAS_FPE
: Type: bool or list of floating-point exception names
: Default: false (unless set by AtlasTool wrapper)
: Effect: Controls enabling floating point exceptions.
: Values:
  - false/0/off/no: disable FPE enabling.
  - true/1/on/yes: enable FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW.
  - Comma-separated list of explicit exception names: FE_INVALID, FE_INEXACT, FE_DIVBYZERO, FE_OVERFLOW, FE_ALL_EXCEPT.
: Source Location: atlas/src/atlas/library/FloatingPointExceptions.cc, atlas/src/atlas/runtime/AtlasTool.cc

ATLAS_SIGNAL_HANDLER
: Type: bool
: Default: false (unless set by AtlasTool wrapper)
: Effect: Enables Atlas signal handlers.
: Values:
  - true/1/on/yes: install Atlas signal handlers.
  - false/0/off/no: do not install.
: Source Location: atlas/src/atlas/library/FloatingPointExceptions.cc, atlas/src/atlas/runtime/AtlasTool.cc

Deprecation Controls
--------------------

ATLAS_DEPRECATION_WARNINGS
: Type: int interpreted as bool (atoi)
: Default: 0
: Effect: Emit deprecation warnings when deprecated factories/builders are used.
: Values:
  - 0: disabled.
  - Non-zero integer: enabled.
: Source Location: atlas/src/atlas/util/Factory.cc

ATLAS_DEPRECATION_ERRORS
: Type: int interpreted as bool (atoi)
: Default: 0
: Effect: Throw an exception on deprecated factory/builder use.
: Values:
  - 0: disabled.
  - Non-zero integer: enabled.
: Source Location: atlas/src/atlas/util/Factory.cc

Interpolation And Geometry
--------------------------

ATLAS_FAST_BUILD_KDTREES
: Type: bool
: Default: true
: Effect: Pre-reserves KD-tree storage in interpolation point-search code paths.
: Values:
  - true/1/on/yes: reserve ahead for faster build in many cases.
  - false/0/off/no: do not pre-reserve.
: Source Location: atlas/src/atlas/interpolation/method/knn/KNearestNeighboursBase.cc, atlas/src/atlas/interpolation/method/PointIndex2.cc, atlas/src/atlas/interpolation/method/PointIndex3.cc, atlas/src/atlas/interpolation/method/PointSet.cc, atlas/src/atlas/interpolation/method/PointSet.h

ATLAS_BUILD_KDTREE_CENTROID
: Type: numeric flag
: Default: 0
: Effect: In conservative spherical polygon interpolation, selects obsolete centroid-based KD-tree build path.
: Values:
  - 0: use standard KD-tree build path.
  - Non-zero: use centroid-based path (marked obsolete in code).
: Source Location: atlas/src/atlas/interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.cc

ATLAS_COMPAREPOINTXYZ_EPS_FACTOR
: Type: double
: Default: 1e8 (multiplied by machine epsilon)
: Effect: Tolerance factor for PointXYZ comparisons in centroid KD-tree path.
: Values:
  - Any positive floating-point number.
: Source Location: atlas/src/atlas/interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.cc

ATLAS_DELAUNAY_BACKEND
: Type: string
: Default: cgal if compiled with CGAL support, otherwise qhull
: Effect: Selects backend for BuildConvexHull3D triangulation.
: Values:
  - qhull: use QHull backend.
  - cgal: use CGAL backend.
: Source Location: atlas/src/atlas/mesh/actions/BuildConvexHull3D.cc

ATLAS_GMSH_FILTER_EDGE_RATIO
: Type: double
: Default: 0.0
: Effect: Gmsh mesh output edge-filter parameter.
: Values:
  - 0.0: no edge-ratio filtering.
  - > 0.0: enable filtering with chosen threshold.
: Source Location: atlas/src/atlas/output/detail/GmshIO.cc

Linear Algebra Backend Selection
--------------------------------

ATLAS_LINALG_FFT_BACKEND
: Type: string
: Default: empty
: Effect: Preferred FFT backend name for Atlas linalg selection.
: Values:
  - Empty: Atlas/eckit default selection.
  - Backend name string (must match a backend available in your build/runtime setup).
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_LINALG_SPARSE_BACKEND
: Type: string
: Default: empty
: Effect: Preferred sparse backend name for Atlas linalg selection.
: Values:
  - Empty: default selection.
  - Backend name string.
: Source Location: atlas/src/atlas/library/Library.cc

ATLAS_LINALG_DENSE_BACKEND
: Type: string
: Default: empty
: Effect: Preferred dense backend name for Atlas linalg selection.
: Values:
  - Empty: default selection.
  - Backend name string.
: Source Location: atlas/src/atlas/library/Library.cc

Debug Selection Variables
-------------------------

ATLAS_DEBUG_GLOBAL_INDEX
: Type: list of gidx_t
: Default: {-1} for accessor, empty for membership checks
: Effect: Selects specific global indices in debug helpers.
: Values:
  - List of integer global indices (typically comma-separated).
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_DEBUG_NODE_GLOBAL_INDEX
: Type: list of gidx_t
: Default: {-1} for accessor, empty for membership checks
: Effect: Selects node global indices in debug helpers.
: Values:
  - List of integer node global indices.
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_DEBUG_EDGE_GLOBAL_INDEX
: Type: list of gidx_t
: Default: {-1} for accessor, empty for membership checks
: Effect: Selects edge global indices in debug helpers.
: Values:
  - List of integer edge global indices.
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_DEBUG_CELL_GLOBAL_INDEX
: Type: list of gidx_t
: Default: {-1} for accessor, empty for membership checks
: Effect: Selects cell global indices in debug helpers.
: Values:
  - List of integer cell global indices.
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_DEBUG_NODE_UID
: Type: list of gidx_t
: Default: {-1} for accessor, empty for membership checks
: Effect: Selects node UIDs in debug helpers.
: Values:
  - List of integer node UIDs.
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_DEBUG_CELL_UID
: Type: list of gidx_t
: Default: empty
: Effect: Selects cell UIDs in debug helpers.
: Values:
  - List of integer cell UIDs.
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_DEBUG_MPI_RANK
: Type: list of long
: Default: {-1} for accessor, empty for membership checks
: Effect: Selects MPI ranks in debug helpers.
: Values:
  - List of integer MPI ranks.
: Source Location: atlas/src/atlas/library/detail/Debug.h

ATLAS_GLOBAL_INDEX
: Type: list of gidx_t
: Default: empty
: Effect: Generic global-index selection in debug helpers.
: Values:
  - List of integer global indices.
: Source Location: atlas/src/atlas/library/detail/Debug.h

Test-Only Variables (atlas/src/tests)
-------------------------------------

ATLAS_MAX_FAILED_EXPECTS
: Type: long
: Default: 100
: Effect: Maximum EXPECT failures allowed before tests abort early.
: Values:
  - Any non-negative integer.
: Source Location: atlas/src/tests/AtlasTestEnvironment.h

ATLAS_MPI_BARRIER_TIMEOUT
: Type: double (seconds)
: Default: 3.0
: Effect: Test helper timeout for MPI barrier deadlock detection.
: Values:
  - Positive floating-point timeout in seconds.
: Source Location: atlas/src/tests/AtlasTestEnvironment.h
