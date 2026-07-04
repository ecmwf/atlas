# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
  TARGET      atlas-main
  OUTPUT_NAME atlas
  SOURCES     atlas.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-meshgen
  SOURCES     atlas-meshgen.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_executable(
  TARGET      atlas-grids
  SOURCES     atlas-grids.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_GRID )

ecbuild_add_executable(
  TARGET      atlas-grid-points
  SOURCES     atlas-grid-points.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-gaussian-latitudes
  SOURCES     atlas-gaussian-latitudes.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_GRID )

ecbuild_add_executable(
  TARGET      atlas-benchmark-relayout
  SOURCES     atlas-benchmark-relayout.cc
  LIBS        atlas )

if( HAVE_TESTS )
  function( atlas_add_benchmark_relayout_smoke test_name )
    ecbuild_add_test(
      TARGET      ${test_name}
      COMMAND     atlas-benchmark-relayout
      ARGS        ${ARGN}
      ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
    )
  endfunction()

  # Small benchmark settings keep these as smoke tests while still exercising code paths.
  set( ATLAS_RELAYOUT_SMOKE_COMMON_ARGS
       --npts=64
       --nproma-other=16
       --iterations=1
       --warmup=0
       --format=json )

  # Rank-4 contiguous paths and nproma dispatch variants.
  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank4_static
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=32
    --loop-order=nproma_innermost
    --nproma-dispatch=static )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank4_static_outer
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=32
    --loop-order=nproma_outermost
    --nproma-dispatch=static )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank4_dynamic
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=32
    --loop-order=nproma_innermost
    --nproma-dispatch=dynamic )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank4_static_default_switch
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=30
    --loop-order=nproma_innermost
    --nproma-dispatch=static )

  # Rank-3 path (nlev > 0, nvar == 0).
  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank3_static
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=4 --nvar=0 --nproma=32
    --loop-order=nproma_innermost
    --nproma-dispatch=static )

  # Rank-2 path (nlev == 0, nvar == 0) with memcpy toggle.
  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank2_memcpy_on
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=0 --nvar=0 --nproma=32
    --blocked-nonblocked-use-memcpy=true )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_rank2_memcpy_off
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=0 --nvar=0 --nproma=32
    --blocked-nonblocked-use-memcpy=false )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_mdspan
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=32
    --implementation=mdspan )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_raw_pointers
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=32
    --implementation=raw_pointers )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_precision_float
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --precision=float
    --nlev=3 --nvar=2 --nproma=32
    --implementation=mdspan )

  atlas_add_benchmark_relayout_smoke( atlas_test_app_benchmark_relayout_blocked_to_blocked_memcpy_off
    ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
    --nlev=3 --nvar=2 --nproma=32
    --blocked-to-blocked-use-memcpy=false )

  ecbuild_add_test(
    TARGET      atlas_test_app_benchmark_relayout_on_device
    COMMAND     atlas-benchmark-relayout
    ARGS        ${ATLAS_RELAYOUT_SMOKE_COMMON_ARGS}
                --nlev=3 --nvar=2 --nproma=32
                --on-device=true
    CONDITION   HAVE_GPU
    ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_RUN_NGPUS=1
  )
  if( TEST atlas_test_app_benchmark_relayout_on_device )
    set_tests_properties( atlas_test_app_benchmark_relayout_on_device PROPERTIES LABELS "gpu" )
  endif()
endif()

add_subdirectory(atlas-interpolations)
