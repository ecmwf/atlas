# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include( ecbuild_add_fortran_flags )

ecbuild_error("Compile flags not tested")

# WARNING!!
# Following has not been thoroughly tested. Take these flags with grain of salt

# Without following line, compile flags are appended to link flags
set( CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( OMP_Fortran_FOUND )
  set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=omp" )
else( )
  set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=noomp" )
endif( )

ecbuild_add_fortran_flags( "-qfree=F90" )
ecbuild_add_fortran_flags( "-qsuffix=cpp=F90" )
ecbuild_add_fortran_flags( "-qextname" )
ecbuild_add_fortran_flags( "-q64=largetype" )
ecbuild_add_fortran_flags( "-qarch=pwr5" )
ecbuild_add_fortran_flags( "-qsource,list" )
ecbuild_add_fortran_flags( "-qsaveopt" )
ecbuild_add_fortran_flags( "-NS32648" )

####################################################################
# RELEASE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O3 -qstrict" BUILD RELEASE)

####################################################################
# DEBUG FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O0 -g" BUILD DEBUG)

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

ecbuild_add_fortran_flags( "-O3 -qstrict" BUILD BIT)

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -b64 -bbigtoc -bmaxstack:0x800000000 -bloadmap:map -bmap:map")

####################################################################

# Meaning of flags
# ----------------
# todo
