include_directories(${GRIDTOOLS_STORAGE_INCLUDE_DIRS})

ecbuild_add_executable(
  TARGET      atlas-main
  OUTPUT_NAME atlas 
  SOURCES     atlas.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-meshgen
  SOURCES     atlas-meshgen.cc
  LIBS        atlas )

ecbuild_add_executable(
  TARGET      atlas-loadbalance
  SOURCES     atlas-loadbalance.cc
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
  LIBS        atlas 
)

ecbuild_add_executable(
  TARGET      atlas-numerics-nabla
  SOURCES     atlas-numerics-nabla.F90
  LIBS        atlas_f
  CONDITION   HAVE_FORTRAN )

