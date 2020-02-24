# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

foreach(test
          test_bounding_box
          test_projection_LAEA
          test_rotation )

    ecbuild_add_test( TARGET atlas_${test} SOURCES ${test}.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )

endforeach()

if(atlas_HAVE_PROJ)
    ecbuild_add_test( TARGET atlas_test_proj SOURCES test_proj.cc LIBS atlas ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} )
endif()

