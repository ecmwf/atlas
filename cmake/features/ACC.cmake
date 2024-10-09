### OpenACC

if( atlas_HAVE_ATLAS_FIELD AND HAVE_GPU )

  if( DEFINED ATLAS_ENABLE_ACC )
    set( ENABLE_ACC ${ATLAS_ENABLE_ACC} )
  endif()
  if( ENABLE_ACC )
    if( NOT HAVE_FORTRAN )
      enable_language(Fortran)
    endif()
    find_package( OpenACC COMPONENTS Fortran CXX )
  endif()

  ecbuild_add_option( FEATURE ACC
                      DESCRIPTION  "OpenACC capable data structures"
                      CONDITION OpenACC_Fortran_FOUND )
  if( HAVE_ACC )
    set( ACC_LINK_OPTIONS ${OpenACC_Fortran_FLAGS} )
  endif()

else()

  set( HAVE_ACC 0 )
  set( atlas_HAVE_ACC 0 )

endif()
