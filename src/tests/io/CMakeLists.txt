# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


ecbuild_add_test( TARGET atlas_test_io_encoding
  SOURCES   test_io_encoding.cc
  LIBS      atlas
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_executable( TARGET atlas_test_io_record
  SOURCES  test_io_record.cc
  LIBS     atlas
  NOINSTALL
)

foreach( algorithm none bzip2 aec lz4 snappy )
    string( TOUPPER ${algorithm} feature )
    if( eckit_HAVE_${feature} OR algorithm MATCHES "none" )
        ecbuild_add_test( TARGET atlas_test_io_record_COMPRESSION_${algorithm}
            COMMAND atlas_test_io_record
            ARGS --suffix ".${algorithm}"
            ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_IO_COMPRESSION=${algorithm}
        )
    endif()
endforeach()

