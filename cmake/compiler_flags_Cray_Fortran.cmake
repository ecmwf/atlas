# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include( ecbuild_add_fortran_flags )

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( NOT ECBUILD_Fortran_FLAGS )
  ecbuild_add_fortran_flags( "-emf -rmoid" )
  ecbuild_add_fortran_flags( "-lhugetlbfs" )
endif()

####################################################################
# RELEASE FLAGS
####################################################################

if( NOT ECBUILD_Fortran_FLAGS )
  ecbuild_add_fortran_flags( "-O3"       BUILD RELEASE )
  ecbuild_add_fortran_flags( "-hfp3"     BUILD RELEASE )
  ecbuild_add_fortran_flags( "-hscalar3" BUILD RELEASE )
  ecbuild_add_fortran_flags( "-hvector3" BUILD RELEASE )
endif()

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

if( NOT ECBUILD_Fortran_FLAGS )
  ecbuild_add_fortran_flags( "-O2"                     BUILD BIT )
  ecbuild_add_fortran_flags( "-hflex_mp=conservative"  BUILD BIT )
  ecbuild_add_fortran_flags( "-hadd_paren"             BUILD BIT )
  ecbuild_add_fortran_flags( "-hfp1"                   BUILD BIT )
endif()

####################################################################
# DEBUG FLAGS
####################################################################

if( NOT ECBUILD_Fortran_FLAGS )
  ecbuild_add_fortran_flags( "-O0"       BUILD DEBUG )
  ecbuild_add_fortran_flags( "-Gfast"    BUILD DEBUG )
  ecbuild_add_fortran_flags( "-Ktrap=fp" BUILD DEBUG )
endif()

####################################################################
# LINK FLAGS
####################################################################

if( NOT ECBUILD_Fortran_LINK_FLAGS )
  set( CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -Wl,-Map,loadmap" )
endif()

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
