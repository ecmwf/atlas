# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

add_subdirectory(test_find_pluto)
add_subdirectory(benchmark)
add_subdirectory(sandbox)

ecbuild_add_test( TARGET pluto_test_pluto_f SOURCES test_pluto_f.F90 LIBS pluto_f )

ecbuild_add_test( TARGET pluto_test_memory_pool SOURCES pluto_test_memory_pool.cc LIBS pluto )
