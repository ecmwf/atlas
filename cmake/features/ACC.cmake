### OpenACC

if( atlas_HAVE_ATLAS_FIELD )

set( ATLAS_ACC_CAPABLE FALSE )
if( HAVE_GPU )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" )
    set( ATLAS_ACC_CAPABLE TRUE )
  else()
    find_package(OpenACC COMPONENTS C Fortran)
    if(OpenACC_Fortran_FOUND AND OpenACC_C_FOUND)
        set( ATLAS_ACC_CAPABLE TRUE )
    endif()
  endif()
endif()

ecbuild_add_option( FEATURE ACC
                    DESCRIPTION  "OpenACC capable data structures"
                    CONDITION ATLAS_ACC_CAPABLE )

if( atlas_HAVE_ACC )
  if( CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC" )
      #set( ACC_Fortran_FLAGS -acc -ta=tesla,nordc )
    set( ACC_Fortran_FLAGS "-acc=gpu;-gpu=gvmode,lineinfo,fastmath,rdc" )
    set( ACC_C_FLAGS ${ACC_Fortran_FLAGS} )
    find_program( ACC_C_COMPILER NAMES pgcc HINTS ${PGI_DIR} ${NVPHC_DIR} ENV PGI_DIR NVHPC_DIR PATH_SUFFIXES bin )
    if( NOT ACC_C_COMPILER )
      ecbuild_error( "Could not find OpenACC capable C compiler" )
    endif()
  else()
    set( ACC_Fortran_FLAGS ${OpenACC_Fortran_FLAGS} )
    set( ACC_C_FLAGS       ${OpenACC_C_FLAGS} )
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

