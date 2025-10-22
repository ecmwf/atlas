cmake_minimum_required( VERSION 3.18 )

project( atlas_field_api VERSION 0.0.0 LANGUAGES C Fortran CXX )

find_package( atlas REQUIRED )
find_package( fckit REQUIRED )
find_package( pluto REQUIRED )
find_package( OpenMP COMPONENTS Fortran )
find_package( OpenACC COMPONENTS Fortran )

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
