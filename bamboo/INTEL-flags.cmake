# don't use trans library until it supports build without MPI
SET( ENABLE_TRANS  ON  CACHE BOOL "Enable TRANS" )
SET( ENABLE_MPI    OFF CACHE BOOL "Disable MPI under Intel compilation" )
