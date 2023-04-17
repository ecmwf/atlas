# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

ecbuild_add_executable(
  TARGET      atlas-main
  OUTPUT_NAME atlas
  SOURCES     atlas.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-meshgen
  SOURCES     atlas-meshgen.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_FUNCTIONSPACE )

ecbuild_add_executable(
  TARGET      atlas-grids
  SOURCES     atlas-grids.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_GRID )

ecbuild_add_executable(
  TARGET      atlas-grid-points
  SOURCES     atlas-grid-points.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-gaussian-latitudes
  SOURCES     atlas-gaussian-latitudes.cc
  LIBS        atlas
  CONDITION   atlas_HAVE_ATLAS_GRID )
