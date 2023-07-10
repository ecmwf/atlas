### Fortran ...

ecbuild_add_option( FEATURE FORTRAN
                    DESCRIPTION "Provide Fortran bindings"
                    REQUIRED_PACKAGES "fckit VERSION 0.6.2 COMPONENTS ECKIT"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE )

if( atlas_HAVE_FORTRAN )

  if( fckit_HAVE_ECKIT )

    ecbuild_enable_fortran( REQUIRED MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/module )
    set( Fortran Fortran )

    if( HAVE_TESTS )
      set( HAVE_FCTEST ON )
      set( atlas_HAVE_FCTEST ON )
    else()
      set( HAVE_FCTEST OFF )
      set( atlas_HAVE_FCTEST OFF )
    endif()

  else()

    ecbuild_warn( "In order to compile atlas_f, fckit is required to be compiled with eckit. Turning off fortran." )
    set( HAVE_FORTRAN 0 )
    set( atlas_HAVE_FORTRAN 0 )

  endif()

endif()

ecbuild_find_python( NO_LIBS )
