# Atlas Environment Variables {#atlas_environment_variables}

This file lists Atlas environment variables found in code under atlas/src (and test-only ones under atlas/src/tests).

## Table of Contents

<ul>
  <li><a href="#notes">Notes</a></li><br>
  <li>
    <a href="#core-logging-and-runtime">Core Logging and Runtime</a><br/>
    <a href="#atlas_info">ATLAS_INFO</a>,
    <a href="#atlas_warning">ATLAS_WARNING</a>,
    <a href="#atlas_trace">ATLAS_TRACE</a>,
    <a href="#atlas_trace_memory">ATLAS_TRACE_MEMORY</a>,
    <a href="#atlas_trace_barriers">ATLAS_TRACE_BARRIERS</a>,
    <a href="#atlas_trace_report">ATLAS_TRACE_REPORT</a>,
    <a href="#atlas_log_rank">ATLAS_LOG_RANK</a>,
    <a href="#atlas_log_file">ATLAS_LOG_FILE</a>,
    <a href="#atlas_workdir">ATLAS_WORKDIR</a>,
    <a href="#atlas_finalises_mpi">ATLAS_FINALISES_MPI</a>,
    <a href="#atlas_plugin_path">ATLAS_PLUGIN_PATH</a>,
    <a href="#atlas_data_path">ATLAS_DATA_PATH</a>,
    <a href="#atlas_cache_path">ATLAS_CACHE_PATH</a>
  </li>

  <li>
    <a href="#floating-point-and-signals">Floating Point and Signals</a><br/>
    <a href="#atlas_fpe">ATLAS_FPE</a>,
    <a href="#atlas_signal_handler">ATLAS_SIGNAL_HANDLER</a>
  </li>

  <li>
    <a href="#deprecation-controls">Deprecation Controls</a><br/>
    <a href="#atlas_deprecation_warnings">ATLAS_DEPRECATION_WARNINGS</a>,
    <a href="#atlas_deprecation_errors">ATLAS_DEPRECATION_ERRORS</a>
  </li>

  <li>
    <a href="#interpolation-and-geometry">Interpolation and Geometry</a><br/>
    <a href="#atlas_fast_build_kdtrees">ATLAS_FAST_BUILD_KDTREES</a>,
    <a href="#atlas_delaunay_backend">ATLAS_DELAUNAY_BACKEND</a>,
    <a href="#atlas_gmsh_filter_edge_ratio">ATLAS_GMSH_FILTER_EDGE_RATIO</a>
  </li>

  <li>
    <a href="#linear-algebra-backend-selection">Linear Algebra Backend Selection</a><br/>
    <a href="#atlas_linalg_fft_backend">ATLAS_LINALG_FFT_BACKEND</a>,
    <a href="#atlas_linalg_sparse_backend">ATLAS_LINALG_SPARSE_BACKEND</a>,
    <a href="#atlas_linalg_dense_backend">ATLAS_LINALG_DENSE_BACKEND</a>
  </li>

  <li>
    <a href="#debug-selection-variables">Debug Selection Variables</a><br/>
    <a href="#atlas_debug_global_index">ATLAS_DEBUG_GLOBAL_INDEX</a>,
    <a href="#atlas_debug_node_global_index">ATLAS_DEBUG_NODE_GLOBAL_INDEX</a>,
    <a href="#atlas_debug_edge_global_index">ATLAS_DEBUG_EDGE_GLOBAL_INDEX</a>,
    <a href="#atlas_debug_cell_global_index">ATLAS_DEBUG_CELL_GLOBAL_INDEX</a>,
    <a href="#atlas_debug_node_uid">ATLAS_DEBUG_NODE_UID</a>,
    <a href="#atlas_debug_cell_uid">ATLAS_DEBUG_CELL_UID</a>,
    <a href="#atlas_debug_mpi_rank">ATLAS_DEBUG_MPI_RANK</a>,
    <a href="#atlas_global_index">ATLAS_GLOBAL_INDEX</a>
  </li>

  <li>
    <a href="#test-only-variables-atlassrctests">Test-only Variables (atlas/src/tests)</a><br/>
    <a href="#atlas_max_failed_expects">ATLAS_MAX_FAILED_EXPECTS</a>,
    <a href="#atlas_mpi_barrier_timeout">ATLAS_MPI_BARRIER_TIMEOUT</a>
  </li>
</ul>

## Notes {#notes}

- Boolean values are parsed by eckit translators. In practice use one of: `0`/`1`, `false`/`true`, `off`/`on`, `no`/`yes`.
- List values use eckit vector parsing (typically comma-separated), for example: `10,20,30`.

## Core Logging and Runtime {#core-logging-and-runtime}

### ATLAS_INFO {#atlas_info}

- Type: `bool`
- Default: `true`
- Effect: Enables Atlas info logging channel.
- Values:
  - `true`/`1`/`on`/`yes`: info logs enabled.
  - `false`/`0`/`off`/`no`: info logs disabled.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_WARNING {#atlas_warning}

- Type: `bool`
- Default: `true`
- Effect: Enables Atlas warning logging channel.
- Values:
  - `true`/`1`/`on`/`yes`: warning logs enabled.
  - `false`/`0`/`off`/`no`: warning logs disabled.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_TRACE {#atlas_trace}

- Type: `bool`
- Default: `false`
- Effect: Enables Atlas trace logging channel.
- Values:
  - `true`/`1`/`on`/`yes`: trace logs enabled.
  - `false`/`0`/`off`/`no`: trace logs disabled.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_TRACE_MEMORY {#atlas_trace_memory}

- Type: `bool`
- Default: `false`
- Effect: Enables memory tracing/reporting path in Atlas tracing.
- Values:
  - `true`/`1`/`on`/`yes`: memory tracing enabled.
  - `false`/`0`/`off`/`no`: memory tracing disabled.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_TRACE_BARRIERS {#atlas_trace_barriers}

- Type: `bool`
- Default: `false`
- Effect: Enables barrier tracing.
- Values:
  - `true`/`1`/`on`/`yes`: barrier tracing enabled.
  - `false`/`0`/`off`/`no`: barrier tracing disabled.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_TRACE_REPORT {#atlas_trace_report}

- Type: `bool`
- Default: `false`
- Effect: Prints trace report during finalise() when tracing is compiled in.
- Values:
  - `true`/`1`/`on`/`yes`: emit final trace report.
  - `false`/`0`/`off`/`no`: do not emit report.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_LOG_RANK {#atlas_log_rank}

- Type: `int`
- Default: `0`
- Effect: MPI rank that keeps normal console logging in several code paths.
- Values:
  - Any integer rank.
- Source Location: `atlas/src/atlas/library/Library.cc`, `atlas/src/atlas/library/FloatingPointExceptions.cc`, `atlas/src/atlas/runtime/AtlasTool.cc`

### ATLAS_LOG_FILE {#atlas_log_file}

- Type: `bool`
- Default: `false`
- Effect: In AtlasTool-based programs, write rank-specific log files.
- Values:
  - `true`/`1`/`on`/`yes`: write log files (name derived from executable/display name and rank).
  - `false`/`0`/`off`/`no`: keep logging policy without per-rank logfile output.
- Source Location: `atlas/src/atlas/runtime/AtlasTool.cc`

### ATLAS_WORKDIR {#atlas_workdir}

- Type: `path string`
- Default: current working directory
- Effect: AtlasTool logfile directory base.
- Values:
  - Any filesystem path.
- Source Location: `atlas/src/atlas/runtime/AtlasTool.cc`

### ATLAS_FINALISES_MPI {#atlas_finalises_mpi}

- Type: `bool`
- Default: `false`
- Effect: If enabled, `atlas::Library::finalise()` calls `atlas::mpi::finalise()`.
- Values:
  - `true`/`1`/`on`/`yes`: Atlas finalise also finalises MPI.
  - `false`/`0`/`off`/`no`: Atlas does not finalise MPI here.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_PLUGIN_PATH {#atlas_plugin_path}

- Type: `string`
- Default: `empty`
- Effect: Extra plugin search path(s) for eckit plugin manager (if supported by linked eckit).
- Values:
  - `Empty`: no extra plugin path.
  - Path list (colon-separated on Unix/macOS).
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_DATA_PATH {#atlas_data_path}

- Type: `path list string`
- Default: `empty` (Atlas also adds built-in fallback path entries)
- Effect: Additional data search path entries.
- Values:
  - Path list (colon-separated on Unix/macOS).
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_CACHE_PATH {#atlas_cache_path}

- Type: `path string`
- Default: `/tmp/cache`
- Effect: Cache location used by `atlas::Library::cachePath()`.
- Values:
  - Any filesystem path.
- Source Location: `atlas/src/atlas/library/Library.cc`

## Floating Point and Signals {#floating-point-and-signals}

### ATLAS_FPE {#atlas_fpe}

- Type: `bool` or list of floating-point exception names
- Default: `false` (unless set by AtlasTool wrapper)
- Effect: Controls enabling floating point exceptions.
- Values:
  - `false`/`0`/`off`/`no`: disable FPE enabling.
  - `true`/`1`/`on`/`yes`: enable `FE_INVALID`, `FE_DIVBYZERO`, `FE_OVERFLOW`.
  - Comma-separated list of explicit exception names: `FE_INVALID`, `FE_INEXACT`, `FE_DIVBYZERO`, `FE_OVERFLOW`, `FE_ALL_EXCEPT`.
- Source Location: `atlas/src/atlas/library/FloatingPointExceptions.cc`, `atlas/src/atlas/runtime/AtlasTool.cc`

### ATLAS_SIGNAL_HANDLER {#atlas_signal_handler}

- Type: `bool`
- Default: `false` (unless set by AtlasTool wrapper)
- Effect: Enables Atlas signal handlers.
- Values:
  - `true`/`1`/`on`/`yes`: install Atlas signal handlers.
  - `false`/`0`/`off`/`no`: do not install.
- Source Location: `atlas/src/atlas/library/FloatingPointExceptions.cc`, `atlas/src/atlas/runtime/AtlasTool.cc`

## Deprecation Controls {#deprecation-controls}

### ATLAS_DEPRECATION_WARNINGS {#atlas_deprecation_warnings}

- Type: `int` interpreted as `bool` (`atoi`)
- Default: `0`
- Effect: Emit deprecation warnings when deprecated factories/builders are used.
- Values:
  - `0`: disabled.
  - Non-zero integer: enabled.
- Source Location: `atlas/src/atlas/util/Factory.cc`

### ATLAS_DEPRECATION_ERRORS {#atlas_deprecation_errors}

- Type: `int` interpreted as `bool` (`atoi`)
- Default: `0`
- Effect: Throw an exception on deprecated factory/builder use.
- Values:
  - `0`: disabled.
  - Non-zero integer: enabled.
- Source Location: `atlas/src/atlas/util/Factory.cc`

## Interpolation and Geometry {#interpolation-and-geometry}

### ATLAS_FAST_BUILD_KDTREES {#atlas_fast_build_kdtrees}

- Type: `bool`
- Default: `true`
- Effect: Pre-reserves KD-tree storage in interpolation point-search code paths.
- Values:
  - `true`/`1`/`on`/`yes`: reserve ahead for faster build in many cases.
  - `false`/`0`/`off`/`no`: do not pre-reserve.
- Source Location: `atlas/src/atlas/interpolation/method/knn/KNearestNeighboursBase.cc`, `atlas/src/atlas/interpolation/method/PointIndex2.cc`, `atlas/src/atlas/interpolation/method/PointIndex3.cc`, `atlas/src/atlas/interpolation/method/PointSet.cc`, `atlas/src/atlas/interpolation/method/PointSet.h`

### ATLAS_DELAUNAY_BACKEND {#atlas_delaunay_backend}

- Type: `string`
- Default: `cgal` if compiled with CGAL support, otherwise `qhull`
- Effect: Selects backend for `BuildConvexHull3D` triangulation.
- Values:
  - `qhull`: use QHull backend.
  - `cgal`: use CGAL backend.
- Source Location: `atlas/src/atlas/mesh/actions/BuildConvexHull3D.cc`

### ATLAS_GMSH_FILTER_EDGE_RATIO {#atlas_gmsh_filter_edge_ratio}

- Type: `double`
- Default: `0.0`
- Effect: Gmsh mesh output edge-filter parameter.
- Values:
  - `0.0`: no edge-ratio filtering.
  - `> 0.0`: enable filtering with chosen threshold.
- Source Location: `atlas/src/atlas/output/detail/GmshIO.cc`

## Linear Algebra Backend Selection {#linear-algebra-backend-selection}

### ATLAS_LINALG_FFT_BACKEND {#atlas_linalg_fft_backend}

- Type: `string`
- Default: `empty`
- Effect: Preferred FFT backend name for Atlas linalg selection.
- Values:
  - `Empty`: Atlas/eckit default selection.
  - Backend name string (must match a backend available in your build/runtime setup).
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_LINALG_SPARSE_BACKEND {#atlas_linalg_sparse_backend}

- Type: `string`
- Default: `empty`
- Effect: Preferred sparse backend name for Atlas linalg selection.
- Values:
  - `Empty`: default selection.
  - Backend name string.
- Source Location: `atlas/src/atlas/library/Library.cc`

### ATLAS_LINALG_DENSE_BACKEND {#atlas_linalg_dense_backend}

- Type: `string`
- Default: `empty`
- Effect: Preferred dense backend name for Atlas linalg selection.
- Values:
  - `Empty`: default selection.
  - Backend name string.
- Source Location: `atlas/src/atlas/library/Library.cc`

## Debug Selection Variables {#debug-selection-variables}

### ATLAS_DEBUG_GLOBAL_INDEX {#atlas_debug_global_index}

- Type: `list of gidx_t`
- Default: `{-1}` for accessor, `empty` for membership checks
- Effect: Selects specific global indices in debug helpers.
- Values:
  - List of integer global indices (typically comma-separated).
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_DEBUG_NODE_GLOBAL_INDEX {#atlas_debug_node_global_index}

- Type: `list of gidx_t`
- Default: `{-1}` for accessor, `empty` for membership checks
- Effect: Selects node global indices in debug helpers.
- Values:
  - List of integer node global indices.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_DEBUG_EDGE_GLOBAL_INDEX {#atlas_debug_edge_global_index}

- Type: `list of gidx_t`
- Default: `{-1}` for accessor, `empty` for membership checks
- Effect: Selects edge global indices in debug helpers.
- Values:
  - List of integer edge global indices.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_DEBUG_CELL_GLOBAL_INDEX {#atlas_debug_cell_global_index}

- Type: `list of gidx_t`
- Default: `{-1}` for accessor, `empty` for membership checks
- Effect: Selects cell global indices in debug helpers.
- Values:
  - List of integer cell global indices.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_DEBUG_NODE_UID {#atlas_debug_node_uid}

- Type: `list of gidx_t`
- Default: `{-1}` for accessor, `empty` for membership checks
- Effect: Selects node UIDs in debug helpers.
- Values:
  - List of integer node UIDs.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_DEBUG_CELL_UID {#atlas_debug_cell_uid}

- Type: `list of gidx_t`
- Default: `empty`
- Effect: Selects cell UIDs in debug helpers.
- Values:
  - List of integer cell UIDs.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_DEBUG_MPI_RANK {#atlas_debug_mpi_rank}

- Type: `list of long`
- Default: `{-1}` for accessor, `empty` for membership checks
- Effect: Selects MPI ranks in debug helpers.
- Values:
  - List of integer MPI ranks.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

### ATLAS_GLOBAL_INDEX {#atlas_global_index}

- Type: `list of gidx_t`
- Default: `empty`
- Effect: Generic global-index selection in debug helpers.
- Values:
  - List of integer global indices.
- Source Location: `atlas/src/atlas/library/detail/Debug.h`

## Test-only Variables (atlas/src/tests) {#test-only-variables-atlassrctests}

### ATLAS_MAX_FAILED_EXPECTS {#atlas_max_failed_expects}

- Type: `long`
- Default: `100`
- Effect: Maximum EXPECT failures allowed before tests abort early.
- Values:
  - Any non-negative integer.
- Source Location: `atlas/src/tests/AtlasTestEnvironment.h`

### ATLAS_MPI_BARRIER_TIMEOUT {#atlas_mpi_barrier_timeout}

- Type: `double` (seconds)
- Default: `3.0`
- Effect: Test helper timeout for MPI barrier deadlock detection.
- Values:
  - Positive floating-point timeout in seconds.
- Source Location: `atlas/src/tests/AtlasTestEnvironment.h`
