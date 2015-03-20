
if( CMAKE_Fortran_COMPILER_LOADED )

ecbuild_add_executable( TARGET atlas-sandbox-fortran-object
  SOURCES fortran_object.F90
)

endif()
