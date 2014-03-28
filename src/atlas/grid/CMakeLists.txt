list( APPEND atlas_grid_srcs
Point3.h
Point3.cc
PointIndex3.h
PointIndex3.cc
GribRead.h
GribRead.cc
GribWrite.h
GribWrite.cc
TriangleIntersection.h
TriangleIntersection.cc
PointSet.h
PointSet.cc
FloatCompare.h
FloatCompare.cc
Tesselation.h
Tesselation.cc
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

# ecbuild_add_library(
# 		TARGET      eckit_tesselation
# 		SOURCES     ${eckit_tesselation_srcs}
# 		CONDITION   HAVE_TESSELATION AND HAVE_GRID AND ATLAS_FOUND AND GRIB_API_FOUND AND HAVE_EIGEN AND Boost_FOUND
# 		INCLUDES    ${ATLAS_INCLUDE_DIRS} ${GRIB_API_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${CGAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
# 		LIBS        eckit ${ATLAS_LIBRARIES} ${GRIB_API_LIBRARIES} ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} )
