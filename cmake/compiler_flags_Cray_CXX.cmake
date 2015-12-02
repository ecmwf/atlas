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
    ecbuild_add_cxx_flags("-lhugetlbfs")
endif()

####################################################################
# RELEASE FLAGS
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
    ecbuild_add_cxx_flags("-O3 -hfp3 -hscalar3 -hvector3 -hPIC -DNDEBUG" BUILD RELEASE)
endif()

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
    ecbuild_add_cxx_flags("-O2 -hflex_mp=conservative -hadd_paren -hfp1 -DNDEBUG" BUILD BIT)
endif()

####################################################################
# DEBUG FLAGS
####################################################################

if( NOT ECBUILD_CXX_FLAGS )
    ecbuild_add_cxx_flags("-O0 -Gfast -Ktrap=fp" BUILD DEBUG)
endif()

####################################################################
# LINK FLAGS
####################################################################

if( NOT ECBUILD_CXX_LINK_FLAGS )
    set( CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -Wl,-Map,loadmap" )
endif()

# This should really go in a toolchain file
ecbuild_info("You should really use a toolchain file for Cray")
#set( CMAKE_CXX_LINK_EXECUTABLE   "<CMAKE_CXX_COMPILER>  <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <OBJECTS>  -o <TARGET> <LINK_LIBRARIES> -Wl,--as-needed,-lmpichf90_cray,--no-as-needed -Wl,-Bdynamic")

####################################################################

# Meaning of flags
# ----------------
# -hfp3     : Special optimisation for floating points
# -Ktrap=fp : Abort on NaN
# -R b      : Bounds checking
# -hflex_mp=conservative -hfp1 : Obtain bit-reproducible results
# -hflex_mp=intolerant -hfp2   : Obtain bit-reproducible results (also)
# -hadd_paren : encourage left to right fp evaluation
# -hscalarN , -hvectorN : optimisation for scalar and vectorisation
# -homp/-hnoomp : Enable/Disable OpenMP
# -rmoi : create compiler listing
