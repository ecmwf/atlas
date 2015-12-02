# don't use trans library until it supports build without MPI
SET( ENABLE_MPI    OFF CACHE BOOL "Disable MPI" )
SET( ENABLE_TRANS  ON  CACHE BOOL "Enable TRANS" )
set( ENABLE_ATLAS_TEST_GRIDSPEC ON CACHE BOOL "Enable atlas_test_gridspec")
