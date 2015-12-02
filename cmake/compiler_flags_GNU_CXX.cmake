# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include( ecbuild_add_cxx_flags )

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
  # nothing to add
endif()

####################################################################
# RELEASE FLAGS
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
    ecbuild_add_cxx_flags( "-O3 -DNDEBUG"  BUILD RELEASE )
endif()

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
    ecbuild_add_cxx_flags("-O2 -DNDEBUG" BUILD BIT)
endif()

####################################################################
# DEBUG FLAGS
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
    ecbuild_add_cxx_flags("-O0 -g -ftrapv" BUILD DEBUG)
endif()

####################################################################
# LINK FLAGS
####################################################################

if( NOT ECBUILD_CXX_LINK_FLAGS )
    # nothing to add
endif()

####################################################################

# Meaning of flags
# ----------------
# todo

