# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( atlas_HAVE_ATLAS_FUNCTIONSPACE )
if( HAVE_FCTEST )
    foreach(test
            fctest_griddistribution
            fctest_grids
            fctest_regional_grids
            fctest_stretchedrotatedgaussiangrid
            fctest_state )

        add_fctest( TARGET atlas_${test} SOURCES ${test}.F90 LINKER_LANGUAGE Fortran LIBS atlas_f ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )

    endforeach()
    add_fctest( TARGET atlas_fctest_unstructuredgrid SOURCES fctest_unstructuredgrid.F90 LINKER_LANGUAGE Fortran LIBS atlas_f CONDITION HAVE_TESSELATION )
endif()
endif()

foreach(test
        test_domain
        test_grid_iterator
        test_stretchedrotatedgaussian
        test_grid_cropping
        test_vertical
        test_spacing
        test_largegrid
        test_grid_hash
        )
    ecbuild_add_test( TARGET atlas_${test} SOURCES ${test}.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )
endforeach()

if( atlas_HAVE_ATLAS_FUNCTIONSPACE )
foreach(test
        test_field
        test_grids
        test_state
        test_cubedsphere
        )
    ecbuild_add_test( TARGET atlas_${test} SOURCES ${test}.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )
endforeach()
endif()

set( _WITH_MPI )
if( eckit_HAVE_MPI )
    set( _WITH_MPI MPI 4 )
endif()

ecbuild_add_test( TARGET  atlas_test_distribution_regular_bands
  ${_WITH_MPI}
  SOURCES test_distribution_regular_bands.cc
  LIBS atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
  CONDITION atlas_HAVE_ATLAS_FUNCTIONSPACE
)


file( GLOB grids ${PROJECT_SOURCE_DIR}/doc/example-grids/*.yml )
if( NOT HAVE_PROJ )
    ecbuild_list_exclude_pattern(
        LIST  grids
        REGEX regional_lambert_azimuthal_equal_area_[3|4].yml )
    ecbuild_list_exclude_pattern(
        LIST  grids
        REGEX regional_polar_stereographic_1.yml )
endif()

foreach( grid ${grids} )
    get_filename_component( grid_name ${grid} NAME_WE )
    set( test_name atlas_test_grid_${grid_name} )
    ecbuild_add_test( TARGET ${test_name}
        COMMAND atlas-grids
        ARGS    ${grid} --check --check-uid --check-boundingbox
        ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
    )
endforeach()
