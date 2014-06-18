
if( CMAKE_Fortran_COMPILER_LOADED )

add_custom_target( atlas_sandbox_fortran_modinc_includes SOURCES
       mod1.h
       mod1.f
       mod2.h
       mod2.f
)

ecbuild_add_library( TARGET atlas_sandbox_fortran_modinc
  CONDITION ON
  SOURCES sb_modinc.F90
)

endif()
