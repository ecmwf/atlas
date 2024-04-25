### OpenACC

if( atlas_HAVE_ATLAS_FIELD )

  set( ATLAS_FIND_OPENACC ${HAVE_CUDA} )
  if (DEFINED ATLAS_ENABLE_ACC)
    if(NOT ATLAS_ENABLE_ACC)
      set( ATLAS_FIND_OPENACC OFF )
    endif()
  elseif (DEFINED ENABLE_ACC)
    if (NOT ENABLE_ACC)
      set( ATLAS_FIND_OPENACC OFF )
    endif()
  endif()

  if (ATLAS_FIND_OPENACC)
    enable_language(C)  # So that find_package(OpenACC) looks C components
    find_package(OpenACC)
  endif()

  ecbuild_add_option( FEATURE ACC
                      DESCRIPTION  "OpenACC capable data structures"
                      CONDITION HAVE_CUDA AND OpenACC_C_FOUND )

else()
  set( HAVE_ACC 0 )
  set( atlas_HAVE_ACC 0 )
endif()

