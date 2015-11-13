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

#ecbuild_add_fortran_flags( "-std=f2008" ) # fortran loc function is non-standard
ecbuild_add_fortran_flags( "-ffree-line-length-none" )

####################################################################
# RELEASE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O3"                BUILD RELEASE )
ecbuild_add_fortran_flags( "-funroll-all-loops" BUILD RELEASE )
ecbuild_add_fortran_flags( "-finline-functions" BUILD RELEASE )
#ecbuild_add_fortran_flags( "-fstack-arrays"     BUILD RELEASE )

####################################################################
# DEBUG FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O0 -g"         BUILD DEBUG )
ecbuild_add_fortran_flags( "-fcheck=bounds" BUILD DEBUG )
ecbuild_add_fortran_flags( "-fbacktrace"    BUILD DEBUG )
ecbuild_add_fortran_flags( "-finit-real=snan" BUILD DEBUG )
ecbuild_add_fortran_flags( "-ffpe-trap=invalid,zero,overflow" BUILD DEBUG )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O2 -g"             BUILD BIT )
ecbuild_add_fortran_flags( "-funroll-all-loops" BUILD BIT )
ecbuild_add_fortran_flags( "-finline-functions" BUILD BIT )
#ecbuild_add_fortran_flags( "-fstack-arrays"     BUILD BIT )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################

# Meaning of flags
# ----------------
# -fstack-arrays     : Allocate automatic arrays on the stack (needs large stacksize!!!)
# -funroll-all-loops : Unroll all loops
# -fcheck=bounds     : Bounds checking

