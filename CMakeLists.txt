
cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( mpdata_shallow_water Fortran )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "MPDATA Shallow Water configuration" )

set(OMP ON CACHE BOOL "Enable/Disable OpenMP" )
set(WARNINGS OFF CACHE BOOL "Enable/Disable all warnings" )
set(GRIB ON CACHE BOOL "Enable/Disable Grib" )


get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if( NOT Fortran_COMPILER_NAME MATCHES "mp.*" AND NOT Fortran_COMPILER_NAME MATCHES "ftn" )
  message( FATAL_ERROR "MPI Compiler needed.\n"
    "Please remove this build, export FC to a MPI compiler such as mpif90 or ftn on Cray and rerun." )
endif( )

# make sure that the default build type is a RELEASE
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug BitReproducible Release." FORCE )
endif( NOT CMAKE_BUILD_TYPE )
string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE )

if( NOT CMAKE_BUILD_TYPE MATCHES "DEBUG" )
  add_definitions( -DNDEBUG )
endif( )

# FFLAGS depend on the compiler
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")
  endif( )
  if( ${WARNINGS} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -fstack-arrays -funroll-all-loops -finline-functions" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace" )
  set( CMAKE_Fortran_FLAGS_BIT     "${CMAKE_Fortran_FLAGS} -O2 -fstack-arrays -funroll-all-loops -finline-functions" )

  # -fstack-arrays     : Allocate automatic arrays on the stack
  # -funroll-all-loops : Unroll all loops
  # -fbounds-check     : Bounds checking
  
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp-stubs")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -ip -unroll -inline -vec-report0" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -g -check bounds -traceback -warn all" )
  set( CMAKE_Fortran_FLAGS_BIT     "${CMAKE_Fortran_FLAGS} -O2 -ip -ipo -unroll -inline -vec-report0" )

elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfree=F90 -qsuffix=cpp=F90 -qextname -qsmp=omp -q64=largetype -qarch=pwr5 -g -qsource,list -qsaveopt -NS32648" )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -qstrict" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0" )
  set( CMAKE_Fortran_FLAGS_BIT     "${CMAKE_Fortran_FLAGS} -O3 -qstrict" )
  set( CMAKE_Fortran_LINK_FLAGS    "-b64 -bbigtoc -bmaxstack:0x800000000 -qsmp=omp")
  set( CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")

elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -emf -rmoi")
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -homp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoomp")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -hfp3 -hmpi1 -hscalar3 -hvector3" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -Gfast -Ktrap=fp -R b" )
  set( CMAKE_Fortran_FLAGS_BIT     "${CMAKE_Fortran_FLAGS} -O2 -hflex_mp=conservative -hadd_paren -hfp1 -hmpi1" )
  # -hfp3     : Special optimisation for floating points
  # -Ktrap=fp : Abort on NaN
  # -R b      : Bounds checking
  # -hflex_mp=conservative -hfp1 : Obtain bit-reproducible results
  # -hflex_mp=intolerant -hfp2   : Obtain bit-reproducible results (also)
  # -hadd_paren : encourage left to right fp evaluation
  # -hscalarN , -hvectorN : optimisation for scalar and vectorisation
  # -homp/-hnoomp : Enable/Disable OpenMP
  # -rmoi : create compiler listing
  set( GRIB_API_DIR "/lus/scratch/ecmwf/esm/grib_api/1.11.0/cray/82" )
endif( )

if( ${GRIB} )
  if( GRIB_API_DIR )
    message( STATUS "Grib enabled : ${GRIB_API_DIR}" )
    add_definitions( -DHAVE_GRIB )
    include_directories(${GRIB_API_DIR}/include)
    link_directories(${GRIB_API_DIR}/lib)
   list( APPEND TPL grib_api_f90 grib_api )
  else( )
    message( STATUS "Grib not found : To enable, run \"cmake . -DGRIB_API_DIR=path/to/grib\"" )
    set( GRIB OFF )
  endif( )
endif( )
include_directories( src )

# Datastructure library
add_library( datastructure
  src/common/parallel_module.F90
  src/common/common_module.F90
  src/datastructure/elements_module.F90
  src/datastructure/lagrangep0_module.F90
  src/datastructure/lagrangep1_module.F90
  src/datastructure/grid_module.F90
  src/datastructure/split_globe_module.F90
  src/datastructure/datastruct_module.F90
  src/io/joanna_module.F90
  src/io/gmsh_module.F90
  src/io/grib_module.F90
)
target_link_libraries( datastructure ${TPL} )


# Shallow water executable
file( COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
add_executable( shallow_water src/shallow_water/shallow_water_module.F90 src/shallow_water/shallow_water.F90 )
target_link_libraries( shallow_water datastructure ${TPL} )

# Test for parallelisation
add_executable( test_sync src/tests/test_sync.F90 )
target_link_libraries( test_sync datastructure )

get_directory_property( defs COMPILE_DEFINITIONS )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "Build   : ${CMAKE_BUILD_TYPE}" )
message( STATUS "Fortran : ${Fortran_COMPILER_NAME} -- ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}" )
message( STATUS "FCFLAGS : ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" )
message( STATUS "LDFLAGS : ${CMAKE_Fortran_LINK_FLAGS}" )
message( STATUS "Defines : ${defs}")
message( STATUS "OpenMP  : ${OMP}")
message( STATUS "Grib    : ${GRIB}")
message( STATUS "--------------------------------------------------------------------" )

