
if( HAVE_FORTRAN )

ecbuild_add_executable( TARGET atlas_example_fortran
  LINKER_LANGUAGE Fortran
  SOURCES example_fortran.F90
  LIBS atlas_f
)

endif()
