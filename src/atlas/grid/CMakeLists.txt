list( APPEND atlas_grid_srcs
Field.cc
Field.h
Gaussian.cc
Gaussian.h
GribRead.cc
GribRead.h
GribWrite.cc
GribWrite.h
Grid.cc
Grid.h
MeshCache.cc
MeshCache.h
LatLon.cc
LatLon.h
Latitudes.cc
Latitudes.h
PointIndex3.cc
PointIndex3.h
PointSet.cc
PointSet.h
ReducedGaussianGrid.h
ReducedGaussianGrid.cc
RegularLatLonGrid.h
RegularLatLonGrid.cc
RegularGaussianGrid.h
RegularGaussianGrid.cc
Tesselation.cc
Tesselation.h
TriangleIntersection.cc
TriangleIntersection.h
Unstructured.cc
Unstructured.h
)

ecbuild_add_library(TARGET atlas_grid
					AUTO_VERSION
					INSTALL_HEADERS ALL
					CONDITION HAVE_GRID AND HAVE_TESSELATION
					SOURCES
						${atlas_grid_srcs}
					HEADER_DESTINATION
						${INSTALL_INCLUDE_DIR}/atlas/grid
					INCLUDES
						"${CGAL_INCLUDE_DIRS}"
					LIBS
						atlas "${ECKIT_LIBRARIES}" "${GRIB_API_LIBRARIES}" "${CGAL_LIBRARY}" "${CGAL_3RD_PARTY_LIBRARIES}" "${GMP_LIBRARIES}" "${MPFR_LIBRARIES}" )
