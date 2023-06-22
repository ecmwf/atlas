### Eigen

if( atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_option( FEATURE EIGEN
                    DESCRIPTION "Use Eigen linear algebra library"
                    REQUIRED_PACKAGES Eigen3 )

if( HAVE_EIGEN AND NOT TARGET Eigen3::Eigen )
    # This is the case for older Eigen versions (e.g. 3.2.0)
    ecbuild_add_library( TARGET atlas_eigen3 TYPE INTERFACE )
    target_include_directories( atlas_eigen3 INTERFACE ${EIGEN3_INCLUDE_DIRS} )
    add_library( Eigen3::Eigen ALIAS atlas_eigen3 )
endif()

else()
    set( HAVE_EIGEN 0 )
    set( atlas_HAVE_EIGEN 0 )
endif()