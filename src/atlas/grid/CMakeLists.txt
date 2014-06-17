list( APPEND atlas_grid_srcs
Gaussian.cc
Gaussian.h
Grid.cc
Grid.h
GridSpec.cc
GridSpec.h
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
ReducedLatLonGrid.h
ReducedLatLonGrid.cc
RegularLatLonGrid.h
RegularLatLonGrid.cc
RotatedLatLonGrid.h
RotatedLatLonGrid.cc
RegularGaussianGrid.h
RegularGaussianGrid.cc
Tesselation.cc
Tesselation.h
Unstructured.cc
Unstructured.h
)

if( HAVE_EIGEN )
  list( APPEND atlas_grid_srcs
    TriangleIntersection.cc
    TriangleIntersection.h
  )
endif()

if( HAVE_GRIB )
  list( APPEND atlas_grid_srcs
    FieldSet.cc
    FieldSet.h
    GribRead.cc
    GribRead.h
    GribWrite.cc
    GribWrite.h
    GridBuilder.h
    GridBuilder.cc
    StackGribFile.cc
    StackGribFile.h
  )
endif()


ecbuild_add_library(TARGET atlas_grid
                    AUTO_VERSION
                    INSTALL_HEADERS ALL
                    CONDITION HAVE_GRID AND HAVE_TESSELATION
                    SOURCES ${atlas_grid_srcs}
                    HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/atlas/grid
                    INCLUDES "${CGAL_INCLUDE_DIRS}"
					LIBS atlas "${ECKIT_LIBRARIES}" "${GRIB_API_LIBRARIES}" "${CGAL_LIBRARIES}" )
