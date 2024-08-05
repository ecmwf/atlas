# (C) Copyright 2024 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_test( TARGET hic_test_dummy
  SOURCES   test_hic_dummy.cc
  LIBS      hic
  DEFINITIONS HIC_BACKEND_DUMMY=1
)
if( HAVE_HIP )
	enable_language(HIP)
      set_source_files_properties( test_hic_dummy.cc  PROPERTIES LANGUAGE HIP )
endif()

if( HAVE_HIP OR HAVE_CUDA )
ecbuild_add_test( TARGET hic_test_hic
  SOURCES   test_hic.cc
  LIBS      hic
)
if( HAVE_HIP )
      set_source_files_properties( test_hic.cc PROPERTIES LANGUAGE HIP )
endif()
endif()

add_subdirectory(test_find_hic)
