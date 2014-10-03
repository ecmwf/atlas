
list( APPEND atlas_grid_srcs
Grid.cc
Grid.h
GridSpec.cc
GridSpec.h
GridSpecParams.cc
GridSpecParams.h
Latitudes.cc
Latitudes.h
MeshCache.cc
MeshCache.h
PointIndex3.cc
PointIndex3.h
PointSet.cc
PointSet.h
PolarStereoGraphic.h
PolarStereoGraphic.cc
ReducedLatLon.cc
ReducedLatLon.h
RegularLatLon.cc
RegularLatLon.h
RotatedLatLon.cc
RotatedLatLon.h
Tesselation.cc
Tesselation.h
TriangleIntersection.cc
TriangleIntersection.h
Unstructured.cc
Unstructured.h
)

if( ECKIT_HAVE_GRIB )
list( APPEND atlas_grid_srcs
FieldSet.cc # to be removed once field set doesn't hold a gribhandle anymore
FieldSet.h
Grib.cc # to be abstracted into a plugin
Grib.h

# contains  include "grib_api.h" for grib_get_gaussian_latitudes()
RegularGG.cc
RegularGG.h
ReducedGG.cc
ReducedGG.h
)
endif()

ecbuild_add_library(TARGET atlas_grid
                    AUTO_VERSION
                    INSTALL_HEADERS ALL
                    CONDITION HAVE_GRID AND HAVE_TESSELATION
                    SOURCES ${atlas_grid_srcs}
                    HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/atlas/grid
                    INCLUDES "${CGAL_INCLUDE_DIRS}"
                    LIBS atlas "${ATLAS_EXTRA_LIBRARIES};${CGAL_LIBRARIES}" )
