
if( CMAKE_Fortran_COMPILER_LOADED )

ecbuild_add_library( TARGET atlas_sandbox_fortran_submodule
  CONDITION OFF
  SOURCES sb_module.F90
)

endif()
