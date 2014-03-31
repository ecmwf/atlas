ecbuild_add_executable(
		TARGET      grib_to_gmsh
		SOURCES     grib_to_gmsh.cc
		LIBS        ${ATLAS_LIBRARIES} ${ECKIT_LIBRARIES} )

