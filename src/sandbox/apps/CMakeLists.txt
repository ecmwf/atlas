# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
  TARGET      atlas-loadbalance
  SOURCES     atlas-loadbalance.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-gmsh-extract
  SOURCES     atlas-gmsh-extract.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-benchmark
  SOURCES     atlas-benchmark.cc
  LIBS        atlas ${OMP_CXX}
)

ecbuild_add_executable(
  TARGET      atlas-numerics-nabla
  SOURCES     atlas-numerics-nabla.F90
  LIBS        atlas_f ${OMP_Fortran}
  LINKER_LANGUAGE Fortran
  CONDITION   HAVE_FORTRAN
)