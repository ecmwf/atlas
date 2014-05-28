# (C) Copyright 1996-2014 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

#######################################################################################
# GNU
#######################################################################################
#Gfortran -W -Wall -pedantic-errors -std=2003 -fbounds-check
#-Werror -ftrace=full
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008")
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-openmp")
  endif( )
  if( ENABLE_WARNINGS )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -finline-functions" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O2 -funroll-all-loops -finline-functions" )

  # -fstack-arrays     : Allocate automatic arrays on the stack
  # -funroll-all-loops : Unroll all loops
  # -fcheck=bounds     : Bounds checking
  
#######################################################################################
# Intel
#######################################################################################
#ifort -check all -warn all,nodec,interfaces -gen interfaces -traceback
#-fpe0 -fpstkchk
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp-stubs")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -unroll -inline -vec-report0 -no-heap-arrays" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -check bounds -traceback -warn all -heap-arrays -fpe-all=0 -fpe:0 -check all" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O2 -ip -ipo -unroll -inline -vec-report0 -no-heap-arrays" )
  set( CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -traceback -fp-trap=common" )


#######################################################################################
# XL (IBM)
#######################################################################################

elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  # Without following line, compile flags are appended to link flags
  set( CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=omp" )
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=omp" )
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qsmp=noomp" )
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -qsmp=noomp" )
  endif( )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfree=F90 -qsuffix=cpp=F90 -qextname -q64=largetype -qarch=pwr5 -g -qsource,list -qsaveopt -NS32648" )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O3 -qstrict" )
  set( CMAKE_Fortran_LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS} -b64 -bbigtoc -bmaxstack:0x800000000 -bloadmap:map -bmap:map")

#######################################################################################
# Cray
#######################################################################################
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -emf -rmoid -lhugetlbfs")
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lhugetlbfs")
  if( ${OMP} )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -homp")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -homp")
  else( )
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoomp")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -hnoomp")
  endif( )
  set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -hfp3 -hmpi1 -hscalar3 -hvector3" )
  set( CMAKE_CXX_FLAGS_RELEASE     "-O3 -hfp3 -hmpi1 -hscalar3 -hvector3" )
  set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Gfast -Ktrap=fp -R b" )
  set( CMAKE_CXX_FLAGS_DEBUG       "-O0 -Gfast -Ktrap=fp" )
  set( CMAKE_Fortran_FLAGS_BIT     "-O2 -hflex_mp=conservative -hadd_paren -hfp1 -hmpi1" )
  set( CMAKE_CXX_FLAGS_BIT         "-O2 -hflex_mp=conservative -hadd_paren -hfp1 -hmpi1" )

  set( CMAKE_Fortran_LINK_FLAGS     "-Wl,-Map,loadmap" )

  # -hfp3     : Special optimisation for floating points
  # -Ktrap=fp : Abort on NaN
  # -R b      : Bounds checking
  # -hflex_mp=conservative -hfp1 : Obtain bit-reproducible results
  # -hflex_mp=intolerant -hfp2   : Obtain bit-reproducible results (also)
  # -hadd_paren : encourage left to right fp evaluation
  # -hscalarN , -hvectorN : optimisation for scalar and vectorisation
  # -homp/-hnoomp : Enable/Disable OpenMP
  # -rmoi : create compiler listing

endif( )

