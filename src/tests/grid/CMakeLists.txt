# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( HAVE_FCTEST )
    foreach(test
            fctest_griddistribution
            fctest_state)

        add_fctest( TARGET atlas_${test} SOURCES ${test}.F90 LINKER_LANGUAGE Fortran LIBS atlas_f )

    endforeach()
endif()

foreach(test
        test_domain
        test_field
        test_grid_ptr
        test_grids
        test_grid_cropping
        test_vertical
        test_state
        test_grid_hash)

    ecbuild_add_test( TARGET atlas_${test} SOURCES ${test}.cc LIBS atlas )

endforeach()
