# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include( ecbuild_check_fortran_source )

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

ecbuild_add_fortran_flags( "-emf -rmoid" )
ecbuild_add_fortran_flags( "-lhugetlbfs" )

####################################################################
# RELEASE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O3"       BUILD RELEASE )
ecbuild_add_fortran_flags( "-hfp3"     BUILD RELEASE )
ecbuild_add_fortran_flags( "-hscalar3" BUILD RELEASE )
ecbuild_add_fortran_flags( "-hvector3" BUILD RELEASE )

####################################################################
# DEBUG FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O0"       BUILD DEBUG )
ecbuild_add_fortran_flags( "-Gfast"    BUILD DEBUG )
ecbuild_add_fortran_flags( "-Ktrap=fp" BUILD DEBUG )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O2"                     BUILD BIT )
ecbuild_add_fortran_flags( "-hflex_mp=conservative"  BUILD BIT )
ecbuild_add_fortran_flags( "-hadd_paren"             BUILD BIT )
ecbuild_add_fortran_flags( "-hfp1"                   BUILD BIT )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "-Wl,-Map,loadmap" )

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
