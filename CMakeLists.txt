
cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

# Disable in source builds (must be called before project)
get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)
if(${srcdir} STREQUAL ${bindir})
  message("######################################################")
  message("You are attempting to build in your Source Directory.")
  message("You must run cmake from a build directory.")
  message("######################################################")
  message(FATAL_ERROR "In-source builds are forbidden!")
endif()
project( pantarhei Fortran )

#######################################################################################
# Configurations
#######################################################################################

set( OMP      ON  CACHE BOOL "Enable/Disable OpenMP"       )
set( WARNINGS OFF CACHE BOOL "Enable/Disable all warnings" )
set( GRIB     ON  CACHE BOOL "Enable/Disable Grib"         )
set( MPI      ON  CACHE BOOL "Enable/Disable MPI"          )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "PantaRhei configuration" )

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if( NOT Fortran_COMPILER_NAME MATCHES "mp.*" AND NOT Fortran_COMPILER_NAME MATCHES "ftn" )
  #message( FATAL_ERROR "MPI Compiler needed.\n"
  #  "Please remove this build, export FC to a MPI compiler such as mpif90 or ftn on Cray and rerun." )
  set( MPI OFF )
endif( )

# make sure that the default build type is a RELEASE
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING 
       "Choose the type of build, options are: Debug Bit Release." FORCE )
endif( )
string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE )

if( NOT CMAKE_BUILD_TYPE MATCHES "DEBUG" )
  add_definitions( -DNDEBUG )
endif( )

if( ${MPI} )
  add_definitions( -DHAVE_MPI )
endif( )

if( ${OMP} )
  add_definitions( -DHAVE_OMP )
endif( )

if( ${GRIB} )
  if( NOT GRIB_API_DIR )
    set( GRIB_API_DIR $ENV{GRIB_API_DIR} )
  endif( )
  if( GRIB_API_DIR )
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

#######################################################################################
# GNU
#######################################################################################
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")
  endif( )
  if( ${WARNINGS} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -finline-functions" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O2 -funroll-all-loops -finline-functions" )

  # -fstack-arrays     : Allocate automatic arrays on the stack
  # -funroll-all-loops : Unroll all loops
  # -fbounds-check     : Bounds checking
  
#######################################################################################
# Intel
#######################################################################################
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp-stubs")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -unroll -inline -vec-report0" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -check bounds -traceback -warn all" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O2 -ip -ipo -unroll -inline -vec-report0" )

#######################################################################################
# XL (IBM)
#######################################################################################
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  # Without following line, compile flags are appended to link flags
  set( CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=omp" )
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=omp" )
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=noomp" )
  endif( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfree=F90 -qsuffix=cpp=F90 -qextname -q64=largetype -qarch=pwr5 -g -qsource,list -qsaveopt -NS32648" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O3 -qstrict" )
  set( CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -b64 -bbigtoc -bmaxstack:0x800000000 -bloadmap:map -bmap:map")

#######################################################################################
# Cray
#######################################################################################
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -emf -rmoid")
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -homp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoomp")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -hfp3 -hmpi1 -hscalar3 -hvector3" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Gfast -Ktrap=fp -R b" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O2 -hflex_mp=conservative -hadd_paren -hfp1 -hmpi1" )
  # -hfp3     : Special optimisation for floating points
  # -Ktrap=fp : Abort on NaN
  # -R b      : Bounds checking
  # -hflex_mp=conservative -hfp1 : Obtain bit-reproducible results
  # -hflex_mp=intolerant -hfp2   : Obtain bit-reproducible results (also)
  # -hadd_paren : encourage left to right fp evaluation
  # -hscalarN , -hvectorN : optimisation for scalar and vectorisation
  # -homp/-hnoomp : Enable/Disable OpenMP
  # -rmoi : create compiler listing
endif( )


#######################################################################################
# Datastructure library
#######################################################################################
if( NOT ${MPI} )
  message( STATUS "MPI not found, use stubs instead" )
  set( MPI_STUBS src/common/mpi_stubs.F90 )
endif( )
add_library( datastructure ${MPI_STUBS}
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

#######################################################################################
# Shallow water executable
#######################################################################################
file( COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
add_executable( shallow_water src/shallow_water/shallow_water_module.F90 src/shallow_water/shallow_water.F90 )
target_link_libraries( shallow_water datastructure ${TPL} )

#######################################################################################
# Test for parallelisation
#######################################################################################
add_executable( test_sync src/tests/test_sync.F90 )
target_link_libraries( test_sync datastructure )


#######################################################################################
# Print Summary
#######################################################################################

get_directory_property( defs COMPILE_DEFINITIONS )

message( STATUS "--------------------------------------------------------------------" )
message( STATUS "System  : ${CMAKE_SYSTEM} -- ${CMAKE_SYSTEM_PROCESSOR}" )
message( STATUS "Build   : ${CMAKE_BUILD_TYPE}" )
message( STATUS "Fortran : ${Fortran_COMPILER_NAME} -- ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}" )
message( STATUS "FCFLAGS : ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}" )
message( STATUS "LDFLAGS : ${CMAKE_Fortran_LINK_FLAGS}" )
message( STATUS "Defines : ${defs}")
message( STATUS "MPI     : ${MPI}")
message( STATUS "OpenMP  : ${OMP}")
message( STATUS "Grib    : ${GRIB}")
message( STATUS "--------------------------------------------------------------------" )

