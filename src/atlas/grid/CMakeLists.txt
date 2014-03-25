list( APPEND atlas_grid_srcs
Field.cc
Field.h
Grid.cc
Grid.h
LatLon.cc
LatLon.h
Gaussian.cc
Gaussian.h
Latitudes.cc
Latitudes.h
 )


ecbuild_add_library(TARGET atlas_grid
					INSTALL_HEADERS ALL
					CONDITION HAVE_GRID
					SOURCES
						${atlas_grid_srcs}
					HEADER_DESTINATION
						${INSTALL_INCLUDE_DIR}/atlas/grid
					LIBS
						atlas ${ECKIT_LIBRARIES} )

