cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project( hello_world VERSION 0.0.0 LANGUAGES CXX )

find_package( atlas REQUIRED )

add_executable( hello_world hello_world.cc )
target_link_libraries( hello_world PUBLIC atlas )

enable_testing()
add_test( NAME test_hello_world COMMAND hello_world )


### Print some information
message( STATUS "atlas_VERSION           ${atlas_VERSION}" )
message( STATUS "atlas_VERSION_STR       ${atlas_VERSION_STR}" )
message( STATUS "atlas_DIR               ${atlas_DIR}" )
message( STATUS "atlas_HAVE_OMP          ${atlas_HAVE_OMP}")
message( STATUS "atlas_HAVE_OMP_CXX      ${atlas_HAVE_OMP_CXX}")
message( STATUS "atlas_HAVE_OMP_Fortran  ${atlas_HAVE_OMP_Fortran}")
message( STATUS "atlas_HAVE_ECTRANS      ${atlas_HAVE_ECTRANS}")
message( STATUS "atlas_HAVE_MPI          ${atlas_HAVE_MPI}")
message( STATUS "ATLAS_LIBRARIES         ${ATLAS_LIBRARIES}")
