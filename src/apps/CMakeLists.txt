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


