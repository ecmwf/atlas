# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include( ecbuild_check_cxx_source )

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( HAVE_OMP )
  cmake_add_cxx_flags( "-fopenmp" )
else( )
  cmake_add_cxx_flags( "-fno-openmp" )
endif( )

####################################################################
# RELEASE FLAGS
####################################################################

cmake_add_cxx_flags( "-O3"  BUILD RELEASE )

####################################################################
# DEBUG FLAGS
####################################################################

cmake_add_cxx_flags( "-O0 -g" BUILD DEBUG )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

cmake_add_cxx_flags( "-O2" BUILD BIT )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_CXX_LINK_FLAGS        "" )

####################################################################

# Meaning of flags
# ----------------
# todo
  
