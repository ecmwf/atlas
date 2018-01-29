
# the following line does not work (CMAKE_SOURCE_DIR seems to point to the build dir)
#include(${CMAKE_SOURCE_DIR}/cmake/atlas_compiler_flags.cmake)

#include(../../git/atlas/cmake/atlas_compiler_flags.cmake)


# on openSUSE 11.3 machines, we need to point to a newer compiler
#SET(CMAKE_Fortran_COMPILER /usr/local/apps/gcc/4.8.1/LP64/bin/gfortran CACHE STRING "Fortran compiler")
#SET(CMAKE_C_COMPILER       /usr/local/apps/gcc/4.8.1/LP64/bin/gcc      CACHE STRING "C compiler")
#SET(CMAKE_CXX_COMPILER     /usr/local/apps/gcc/4.8.1/LP64/bin/g++      CACHE STRING "C++ compiler")
#set(CMAKE_Fortran_LINK_FLAGS  "-L/usr/local/apps/gcc/4.8.1/LP64/lib/gcc/x86_64-suse-linux/4.8.1/")
#link_directories("/usr/local/apps/gcc/4.8.1/LP64/lib/gcc/x86_64-suse-linux/4.8.1/")

SET(ENABLE_SANDBOX  OFF  CACHE BOOL "Disable sandbox")
