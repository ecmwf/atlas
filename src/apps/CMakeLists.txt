ecbuild_add_executable(
  TARGET      atlas-main
  OUTPUT_NAME atlas 
  SOURCES     atlas.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      grib_to_gmsh
  SOURCES     grib_to_gmsh.cc
  LIBS        atlas
  CONDITION   HAVE_GRID AND HAVE_GRIB )

ecbuild_add_executable(
  TARGET      gmsh_ls
  SOURCES     gmsh_ls.cc
  LIBS        atlas
  CONDITION   HAVE_GRID )

ecbuild_add_executable(
  TARGET      atlas-meshgen
  SOURCES     atlas-meshgen.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-gmsh-extract
  SOURCES     atlas-gmsh-extract.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-grids
  SOURCES     atlas-grids.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-gaussian-latitudes
  SOURCES     atlas-gaussian-latitudes.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-benchmark
  SOURCES     atlas-benchmark.cc
  LIBS        atlas )
