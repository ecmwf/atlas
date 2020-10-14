### OpenACC

set( ATLAS_ACC_CAPABLE FALSE )
if( HAVE_CUDA )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
    set( ATLAS_ACC_CAPABLE TRUE )
  endif()
endif()

ecbuild_add_option( FEATURE ACC
                    DESCRIPTION  "OpenACC capable data structures"
                    CONDITION ATLAS_ACC_CAPABLE )

if( atlas_HAVE_ACC )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
    set( ACC_Fortran_FLAGS -acc -ta=tesla,nordc )
    set( ACC_C_FLAGS ${ACC_Fortran_FLAGS} )
    find_program( ACC_C_COMPILER NAMES pgcc HINTS ${PGI_DIR} ENV PGI_DIR PATH_SUFFIXES bin )
    if( NOT ACC_C_COMPILER )
      ecbuild_error( "Could not find OpenACC capable C compiler" )
    endif()
  endif()
endif()
