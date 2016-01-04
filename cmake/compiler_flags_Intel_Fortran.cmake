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
  # Nothing to add
endif()

####################################################################
# RELEASE FLAGS
####################################################################

if( NOT ECBUILD_Fortran_FLAGS_RELEASE )
  ecbuild_add_fortran_flags( "-O3"             BUILD RELEASE )
  ecbuild_add_fortran_flags( "-unroll"         BUILD RELEASE )
  ecbuild_add_fortran_flags( "-inline"         BUILD RELEASE )
  ecbuild_add_fortran_flags( "-vec-report0"    BUILD RELEASE )
  ecbuild_add_fortran_flags( "-heap-arrays"    BUILD RELEASE )
endif()

####################################################################
# DEBUG FLAGS
####################################################################

if( NOT ECBUILD_Fortran_FLAGS_DEBUG )
  ecbuild_add_fortran_flags( "-O0 -g"                       BUILD DEBUG )
  ecbuild_add_fortran_flags( "-check bounds"                BUILD DEBUG )
  ecbuild_add_fortran_flags( "-traceback"                   BUILD DEBUG )
  ecbuild_add_fortran_flags( "-warn all"                    BUILD DEBUG )
  ecbuild_add_fortran_flags( "-heap-arrays"                 BUILD DEBUG )
  ecbuild_add_fortran_flags( "-fpe-all=0 -fpe:0 -check all" BUILD DEBUG )
endif()

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

if( NOT ECBUILD_Fortran_FLAGS_BIT )
  ecbuild_add_fortran_flags( "-O2"             BUILD BIT )
  ecbuild_add_fortran_flags( "-unroll"         BUILD BIT )
  ecbuild_add_fortran_flags( "-inline"         BUILD BIT )
  ecbuild_add_fortran_flags( "-vec-report0"    BUILD BIT )
  ecbuild_add_fortran_flags( "-heap-arrays"    BUILD BIT )
endif()

####################################################################
# LINK FLAGS
####################################################################

if( NOT ECBUILD_Fortran_LINK_FLAGS )
  # nothing to add
endif()

####################################################################

# Meaning of flags
# ----------------
# todo

