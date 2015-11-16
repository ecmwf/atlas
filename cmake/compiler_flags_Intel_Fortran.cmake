# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################


####################################################################
# RELEASE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O3"             BUILD RELEASE )
ecbuild_add_fortran_flags( "-unroll"         BUILD RELEASE )
ecbuild_add_fortran_flags( "-inline"         BUILD RELEASE )
ecbuild_add_fortran_flags( "-vec-report0"    BUILD RELEASE )
ecbuild_add_fortran_flags( "-no-heap-arrays" BUILD RELEASE )

####################################################################
# DEBUG FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O0 -g"                       BUILD DEBUG )
ecbuild_add_fortran_flags( "-check bounds"                BUILD DEBUG )
ecbuild_add_fortran_flags( "-traceback"                   BUILD DEBUG )
ecbuild_add_fortran_flags( "-warn all"                    BUILD DEBUG )
ecbuild_add_fortran_flags( "-heap-arrays"                 BUILD DEBUG )
ecbuild_add_fortran_flags( "-fpe-all=0 -fpe:0 -check all" BUILD DEBUG )

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O2"             BUILD BIT )
ecbuild_add_fortran_flags( "-unroll"         BUILD BIT )
ecbuild_add_fortran_flags( "-inline"         BUILD BIT )
ecbuild_add_fortran_flags( "-vec-report0"    BUILD BIT )
ecbuild_add_fortran_flags( "-no-heap-arrays" BUILD BIT )

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "" )

####################################################################

# Meaning of flags
# ----------------
# todo

