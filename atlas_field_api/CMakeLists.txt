cmake_minimum_required( VERSION 3.18 )

project( atlas_field_api VERSION 0.0.0 LANGUAGES C Fortran CXX )

find_package( atlas REQUIRED )
find_package( fckit REQUIRED )
find_package( pluto REQUIRED )
find_package( OpenMP COMPONENTS Fortran )

set (HAVE_ACC 0)
if (atlas_HAVE_ACC)
  find_package( OpenACC COMPONENTS Fortran REQUIRED )
  set (HAVE_ACC 1)
endif()
message(STATUS "HAVE_ACC: ${HAVE_ACC}")

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
