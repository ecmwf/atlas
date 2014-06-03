# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

include( ecbuild_check_fortran_source )

# WARNING!!
# Following has not been thoroughly tested. Take these flags with grain of salt


# Without following line, compile flags are appended to link flags
set( CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")

####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

if( ${HAVE_OMP} )
  cmake_add_fortran_flags( "-qsmp=omp" )
  set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=omp" )
else( )
  cmake_add_fortran_flags( "-qsmp=noomp" )
  set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=noomp" )
endif( )

cmake_add_fortran_flags( "-qfree=F90" )
cmake_add_fortran_flags( "-qsuffix=cpp=F90" )
cmake_add_fortran_flags( "-qextname" )
cmake_add_fortran_flags( "-q64=largetype" )
cmake_add_fortran_flags( "-qarch=pwr5" )
cmake_add_fortran_flags( "-qsource,list" )
cmake_add_fortran_flags( "-qsaveopt" )
cmake_add_fortran_flags( "-NS32648" )

####################################################################
# RELEASE FLAGS
####################################################################

cmake_add_fortran_flags( "-O3 -qstrict" BUILD RELEASE)

####################################################################
# DEBUG FLAGS
####################################################################

cmake_add_fortran_flags( "-O0 -g" BUILD DEBUG)

####################################################################
# BIT REPRODUCIBLE FLAGS
####################################################################

cmake_add_fortran_flags( "-O3 -qstrict" BUILD BIT)

####################################################################
# LINK FLAGS
####################################################################

set( CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -b64 -bbigtoc -bmaxstack:0x800000000 -bloadmap:map -bmap:map")

####################################################################

# Meaning of flags
# ----------------
# todo
