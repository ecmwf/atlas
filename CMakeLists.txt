
cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( mpdata_shallow_water Fortran )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "MPDATA Shallow Water configuration" )

# make sure that the default build type is a RELEASE
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
endif( NOT CMAKE_BUILD_TYPE )
string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE )

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-fstack-arrays -funroll-all-loops -fno-f2c -O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -hfp3" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Ktrap=fp -R b" )
endif( )

if( NOT Fortran_COMPILER_NAME MATCHES "ftn" )
  find_package( mpi )
  include_directories( ${MPI_Fortran_INCLUDE_PATH} )
endif( )

# Datastructure library
add_library( datastructure
  src/common/parallel_module.F90
  src/common/common_module.F90
  src/mesh/elements_module.F90
  src/mesh/lagrangep0_module.F90
  src/mesh/lagrangep1_module.F90
  src/mesh/grid_module.F90
  src/mesh/split_globe_module.F90
  src/mesh/datastruct_module.F90
  src/io/read_joanna_module.F90
  src/io/gmsh_module.F90
  src/io/grib_module.F90
)
target_link_libraries( datastructure ${MPI_Fortran_LIBRARIES} )

# Shallow water executable
add_executable( shallow_water src/shallow_water_module.F90 src/shallow_water.F90 )
target_link_libraries( shallow_water datastructure )

# Test for parallelisation
add_executable( test_sync src/test_sync.F90 )
target_link_libraries( test_sync datastructure )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "Build   : ${CMAKE_BUILD_TYPE}" )
message( STATUS "Fortran : ${Fortran_COMPILER_NAME} (${CMAKE_Fortran_COMPILER_ID})" )
message( STATUS "Flags   : ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" )
message( STATUS "--------------------------------------------------------------------" )

