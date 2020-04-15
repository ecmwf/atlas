cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project( hello_world VERSION 0.0.0 LANGUAGES Fortran )

find_package( atlas REQUIRED )

if( NOT atlas_HAVE_FORTRAN )
  message( FATAL_ERROR "Atlas was not compiled with Fortran enabled" )
endif()


add_library( hw SHARED hello_world_module.F90 )
target_link_libraries( hw PUBLIC atlas_f )

add_executable( hello_world hello_world.F90 )
target_link_libraries( hello_world PUBLIC hw )
set_target_properties( hello_world PROPERTIES LINKER_LANGUAGE Fortran )

enable_testing()
add_test( NAME test_hello_world COMMAND hello_world )
