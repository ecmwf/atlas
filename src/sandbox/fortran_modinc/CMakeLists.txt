
if( CMAKE_Fortran_COMPILER_LOADED )

ecbuild_add_library( TARGET atlas_sandbox_fortran_submodule
  CONDITION ON
  SOURCES sb_modinc.F90
)

endif()
