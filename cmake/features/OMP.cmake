### OMP ...
if( "${CMAKE_VERSION}" VERSION_LESS "3.11" )
  if( ENABLE_OMP )
    ecbuild_warn( "OpenMP only supported with CMake 3.11 onwards" )
  endif()
else()
  find_package( OpenMP COMPONENTS CXX ${Fortran} )
endif()
ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "support for OpenMP shared memory parallelism"
                    CONDITION OpenMP_Fortran_FOUND OR OpenMP_CXX_FOUND )
ecbuild_add_option( FEATURE OMP_Fortran
                    DESCRIPTION "support for Fortran OpenMP shared memory parallelism"
                    CONDITION HAVE_OMP AND OpenMP_Fortran_FOUND )

ecbuild_add_option( FEATURE OMP_CXX
                    DESCRIPTION "support for CXX OpenMP shared memory parallelism"
                    CONDITION HAVE_OMP AND OpenMP_CXX_FOUND )
if( TARGET OpenMP::OpenMP_CXX )
  set( OMP_CXX OpenMP::OpenMP_CXX )
endif()
if( TARGET OpenMP::OpenMP_Fortran )
  set( OMP_Fortran OpenMP::OpenMP_Fortran )
endif()
