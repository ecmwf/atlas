# don't use trans library until it supports build without MPI
set( ENABLE_MPI                 OFF CACHE BOOL "Disable MPI under GNU compilation" )
set( ENABLE_TRANS               ON  CACHE BOOL "Enable TRANS" )
set( ENABLE_ATLAS_TEST_GRIDSPEC ON  CACHE BOOL "Enable atlas_test_gridspec")
set( ENABLE_BOUNDSCHECKING      ON  CACHE BOOL "Enable bounds checking")

set(CMAKE_Fortran_FLAGS "-ffree-line-length-none" CACHE STRING "Fortran compiler flags")
