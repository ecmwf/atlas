# (C) Copyright 2022 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET atlas_io_test_encoding
  SOURCES   test_io_encoding.cc
  LIBS      atlas_io
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_test( TARGET atlas_io_test_stream
  SOURCES   test_io_stream.cc
  LIBS      atlas_io
  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
)

ecbuild_add_executable( TARGET atlas_io_test_record
  SOURCES  test_io_record.cc
  LIBS     atlas_io
  NOINSTALL
)

foreach( algorithm none bzip2 aec lz4 snappy )
    string( TOUPPER ${algorithm} feature )
    if( eckit_HAVE_${feature} OR algorithm MATCHES "none" )
        ecbuild_add_test( TARGET atlas_io_test_record_COMPRESSION_${algorithm}
            COMMAND atlas_io_test_record
            ARGS --suffix ".${algorithm}"
            ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT} ATLAS_IO_COMPRESSION=${algorithm}
        )
    endif()
endforeach()

