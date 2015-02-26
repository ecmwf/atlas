# don't use trans library until it supports build without MPI
SET( ENABLE_TRANS  OFF CACHE BOOL "Disable TRANS library" )
SET( ENABLE_MPI    OFF CACHE BOOL "Disable MPI under Clang compilation" )
