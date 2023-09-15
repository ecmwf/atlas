if( atlas_HAVE_ATLAS_FUNCTIONSPACE )
### tesselation ...

ecbuild_add_option( FEATURE TESSELATION
                    DESCRIPTION "Support for unstructured mesh generation"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE
                    REQUIRED_PACKAGES "Qhull" )
if(HAVE_TESSELATION)
  set(QHULL_LIBRARIES Qhull::qhullcpp Qhull::qhull_r)
  set(atlas_HAVE_QHULL 1)
else()
  set(atlas_HAVE_QHULL 0)
endif()

### NOTE
#
# CGAL is deprecated as TESSELATION backend. Qhull is to be used instead.
# To use CGAL regardless, turn ON CGAL feature (-DENABLE_CGAL=ON)

set(Boost_USE_MULTITHREADED ON )
ecbuild_add_option( FEATURE CGAL
                    DEFAULT OFF
                    DESCRIPTION "Support for unstructured mesh generation"
                    CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE
                    REQUIRED_PACKAGES
                      "CGAL"
                      "Boost VERSION 1.45.0 QUIET" )

if( HAVE_CGAL )
    list( APPEND CGAL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
    if ( TARGET CGAL::CGAL )
      list( APPEND CGAL_LIBRARIES CGAL::CGAL ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
      set( _cgal_target CGAL::CGAL )
      get_target_property(_aliased CGAL::CGAL ALIASED_TARGET)
      if(_aliased)
        set( _cgal_target ${_aliased} )
      endif()
      # Reset INTERFACE_COMPILE_OPTIONS ( see ATLAS-193 )
      get_target_property( CGAL_COMPILE_FLAGS ${_cgal_target} INTERFACE_COMPILE_OPTIONS )
      set_target_properties( ${_cgal_target} PROPERTIES INTERFACE_COMPILE_OPTIONS "" )
    else()
      list( APPEND CGAL_LIBRARIES ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
    endif()
endif()

if( NOT HAVE_CGAL )
    unset( CGAL_LIBRARIES )
    unset( CGAL_INCLUDE_DIRS )
endif()
endif()
