ecbuild_add_executable(
  TARGET      atlas-main
  SOURCES     atlas.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )
set_target_properties( atlas-main
  PROPERTIES
  OUTPUT_NAME atlas )

ecbuild_add_executable(
  TARGET      grib_to_gmsh
  SOURCES     grib_to_gmsh.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES}
  CONDITION   HAVE_GRID AND HAVE_GRIB )

ecbuild_add_executable(
  TARGET      gmsh_ls
  SOURCES     gmsh_ls.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES}
  CONDITION   HAVE_GRID )

ecbuild_add_executable(
  TARGET      atlas-meshgen
  SOURCES     atlas-meshgen.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )

ecbuild_add_executable(
  TARGET      atlas-gmsh-extract
  SOURCES     atlas-gmsh-extract.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )

ecbuild_add_executable(
  TARGET      atlas-grids
  SOURCES     atlas-grids.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )

ecbuild_add_executable(
  TARGET      atlas-gaussian-latitudes
  SOURCES     atlas-gaussian-latitudes.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )

ecbuild_add_executable(
  TARGET      atlas-benchmark
  SOURCES     atlas-benchmark.cc
  LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )
