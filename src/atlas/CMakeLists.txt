# (C) Copyright 1996-2017 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

### config headers

ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( library/defines.h.in   library/defines.h   )
configure_file( library/git_sha1.h.in  library/git_sha1.h   )

install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/library/defines.h
  ${CMAKE_CURRENT_BINARY_DIR}/library/config.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas/build )

### sources


list( APPEND atlas_srcs
library/atlas.h
library/atlas.cc
library/config.h
library/config.cc
library/LibAtlas.h
library/LibAtlas.cc
runtime/ErrorHandling.cc
runtime/ErrorHandling.cc
runtime/AtlasTool.h
runtime/AtlasTool.cc
runtime/Log.h
)

list( APPEND atlas_grid_srcs

grid.h
grid/Grid.cc
grid/Grid.h
grid/Distribution.cc
grid/Distribution.h
grid/Projection.cc
grid/Projection.h
grid/Domain.cc
grid/Domain.h
grid/Spacing.cc
grid/Spacing.h
grid/Partitioner.h
grid/Partitioner.cc
grid/Iterator.h

grid/detail/grid/GridBuilder.h
grid/detail/grid/GridBuilder.cc
grid/detail/grid/Grid.h
grid/detail/grid/Grid.cc
grid/detail/grid/Structured.cc
grid/detail/grid/Structured.h
grid/detail/grid/Unstructured.cc
grid/detail/grid/Unstructured.h
grid/detail/grid/Gaussian.h
grid/detail/grid/Gaussian.cc
grid/detail/grid/LonLat.h
grid/detail/grid/LonLat.cc
grid/detail/grid/Regional.h
grid/detail/grid/Regional.cc

grid/detail/domain/Domain.cc
grid/detail/domain/Domain.h
grid/detail/domain/RectangularDomain.h
grid/detail/domain/RectangularDomain.cc
grid/detail/domain/CircularDomain.h
grid/detail/domain/CircularDomain.cc
grid/detail/domain/EmptyDomain.h
grid/detail/domain/EmptyDomain.cc
grid/detail/domain/GlobalDomain.h
grid/detail/domain/GlobalDomain.cc
grid/detail/domain/ZonalBandDomain.h
grid/detail/domain/ZonalBandDomain.cc

grid/detail/partitioner/EqualRegionsPartitioner.cc
grid/detail/partitioner/EqualRegionsPartitioner.h
grid/detail/partitioner/CheckerboardPartitioner.cc
grid/detail/partitioner/CheckerboardPartitioner.h
grid/detail/partitioner/PartitionerFromPrePartitionedMesh.h
grid/detail/partitioner/PrePartitionedBruteForce.cc
grid/detail/partitioner/PrePartitionedBruteForce.h
grid/detail/partitioner/PrePartitionedPolygon.cc
grid/detail/partitioner/PrePartitionedPolygon.h
grid/detail/partitioner/Partitioner.cc
grid/detail/partitioner/Partitioner.h

grid/detail/projection/Projection.cc
grid/detail/projection/Projection.h
grid/detail/projection/Rotation.cc
grid/detail/projection/Rotation.h
grid/detail/projection/LonLatProjection.cc
grid/detail/projection/LonLatProjection.h
grid/detail/projection/LambertProjection.cc
grid/detail/projection/LambertProjection.h
grid/detail/projection/SchmidtProjection.cc
grid/detail/projection/SchmidtProjection.h
grid/detail/projection/MercatorProjection.cc
grid/detail/projection/MercatorProjection.h

grid/detail/spacing/Spacing.cc
grid/detail/spacing/Spacing.h
grid/detail/spacing/CustomSpacing.h
grid/detail/spacing/CustomSpacing.cc
grid/detail/spacing/LinearSpacing.h
grid/detail/spacing/LinearSpacing.cc
grid/detail/spacing/FocusSpacing.h
grid/detail/spacing/FocusSpacing.cc
grid/detail/spacing/GaussianSpacing.h
grid/detail/spacing/GaussianSpacing.cc
grid/detail/spacing/gaussian/Latitudes.cc
grid/detail/spacing/gaussian/Latitudes.h
grid/detail/spacing/gaussian/N.cc
grid/detail/spacing/gaussian/N.h
grid/detail/spacing/gaussian/N16.cc    # TL31
grid/detail/spacing/gaussian/N24.cc    # TL47
grid/detail/spacing/gaussian/N32.cc    # TL63
grid/detail/spacing/gaussian/N48.cc    # TL95
grid/detail/spacing/gaussian/N64.cc    # TL127
grid/detail/spacing/gaussian/N80.cc    # TL159
grid/detail/spacing/gaussian/N96.cc    # TL191
grid/detail/spacing/gaussian/N128.cc   # TL255
grid/detail/spacing/gaussian/N160.cc   # TL319
grid/detail/spacing/gaussian/N200.cc   # TL399
grid/detail/spacing/gaussian/N256.cc   # TL511
grid/detail/spacing/gaussian/N320.cc   # TL639
grid/detail/spacing/gaussian/N400.cc   # TL799
grid/detail/spacing/gaussian/N512.cc   # TL1023
grid/detail/spacing/gaussian/N576.cc   # TL1151
grid/detail/spacing/gaussian/N640.cc   # TL1279
grid/detail/spacing/gaussian/N800.cc   # TL1599
grid/detail/spacing/gaussian/N1024.cc  # TL2047
grid/detail/spacing/gaussian/N1280.cc  # TL2559
grid/detail/spacing/gaussian/N1600.cc  # TL3199
grid/detail/spacing/gaussian/N2000.cc  # TL3999
grid/detail/spacing/gaussian/N4000.cc  # TL7999
grid/detail/spacing/gaussian/N8000.cc  # TL15999

grid/detail/pl/classic_gaussian/N.h
grid/detail/pl/classic_gaussian/N.cc
grid/detail/pl/classic_gaussian/PointsPerLatitude.h
grid/detail/pl/classic_gaussian/PointsPerLatitude.cc
grid/detail/pl/classic_gaussian/N16.cc    # TL31
grid/detail/pl/classic_gaussian/N24.cc    # TL47
grid/detail/pl/classic_gaussian/N32.cc    # TL63
grid/detail/pl/classic_gaussian/N48.cc    # TL95
grid/detail/pl/classic_gaussian/N64.cc    # TL127
grid/detail/pl/classic_gaussian/N80.cc    # TL159
grid/detail/pl/classic_gaussian/N96.cc    # TL191
grid/detail/pl/classic_gaussian/N128.cc   # TL255
grid/detail/pl/classic_gaussian/N160.cc   # TL319
grid/detail/pl/classic_gaussian/N200.cc   # TL399
grid/detail/pl/classic_gaussian/N256.cc   # TL511
grid/detail/pl/classic_gaussian/N320.cc   # TL639
grid/detail/pl/classic_gaussian/N400.cc   # TL799
grid/detail/pl/classic_gaussian/N512.cc   # TL1023
grid/detail/pl/classic_gaussian/N576.cc   # TL1151
grid/detail/pl/classic_gaussian/N640.cc   # TL1279
grid/detail/pl/classic_gaussian/N800.cc   # TL1599
grid/detail/pl/classic_gaussian/N1024.cc  # TL2047
grid/detail/pl/classic_gaussian/N1280.cc  # TL2559
grid/detail/pl/classic_gaussian/N1600.cc  # TL3199
grid/detail/pl/classic_gaussian/N2000.cc  # TL3999
grid/detail/pl/classic_gaussian/N4000.cc  # TL7999
grid/detail/pl/classic_gaussian/N8000.cc  # TL15999
)
if( ATLAS_HAVE_TRANS )
list( APPEND atlas_grid_srcs
  grid/detail/partitioner/TransPartitioner.h
  grid/detail/partitioner/TransPartitioner.cc
)
endif()

list( APPEND atlas_mesh_srcs
mesh/Connectivity.h
mesh/Connectivity.cc
mesh/Elements.h
mesh/Elements.cc
mesh/ElementType.h
mesh/ElementType.cc
mesh/Halo.h
mesh/Halo.cc
mesh/HybridElements.h
mesh/HybridElements.cc
mesh/Mesh.h
mesh/Mesh.cc
mesh/Nodes.h
mesh/Nodes.cc
mesh/actions/ExtendNodesGlobal.h
mesh/actions/ExtendNodesGlobal.cc
mesh/actions/BuildDualMesh.h
mesh/actions/BuildDualMesh.cc
mesh/actions/BuildEdges.h
mesh/actions/BuildEdges.cc
mesh/actions/BuildHalo.h
mesh/actions/BuildHalo.cc
mesh/actions/BuildParallelFields.h
mesh/actions/BuildParallelFields.cc
mesh/actions/BuildPeriodicBoundaries.h
mesh/actions/BuildPeriodicBoundaries.cc
mesh/actions/BuildStatistics.h
mesh/actions/BuildStatistics.cc
mesh/actions/BuildXYZField.h
mesh/actions/BuildXYZField.cc
mesh/actions/BuildCellCentres.h
mesh/actions/BuildCellCentres.cc
mesh/actions/BuildConvexHull3D.h
mesh/actions/BuildConvexHull3D.cc
mesh/actions/GenerateMesh.h
mesh/actions/GenerateMesh.cc
mesh/actions/WriteLoadBalanceReport.h
mesh/actions/WriteLoadBalanceReport.cc
meshgenerator/DelaunayMeshGenerator.h
meshgenerator/DelaunayMeshGenerator.cc
meshgenerator/MeshGenerator.cc
meshgenerator/MeshGenerator.h
meshgenerator/StructuredMeshGenerator.h
meshgenerator/StructuredMeshGenerator.cc
meshgenerator/RegularMeshGenerator.cc
meshgenerator/RegularMeshGenerator.h
mesh/actions/BuildTorusXYZField.h
mesh/actions/BuildTorusXYZField.cc
)

list( APPEND atlas_output_srcs
output/Output.h
output/Output.cc
output/Gmsh.h
output/Gmsh.cc
)

list( APPEND atlas_field_srcs
field/Field.h
field/Field.cc
field/FieldCreator.h
field/FieldCreator.cc
field/FieldCreatorArraySpec.h
field/FieldCreatorArraySpec.cc
field/FieldCreatorIFS.h
field/FieldCreatorIFS.cc
field/FieldSet.h
field/FieldSet.cc
field/State.h
field/State.cc
field/Options.h
)

list( APPEND atlas_functionspace_srcs
functionspace/EdgeColumns.h
functionspace/EdgeColumns.cc
functionspace/FunctionSpace.h
functionspace/FunctionSpace.cc
functionspace/NodeColumns.h
functionspace/NodeColumns.cc
functionspace/NodeColumnsInterface.h
functionspace/NodeColumnsInterface.cc
functionspace/StructuredColumns.h
functionspace/StructuredColumns.cc
functionspace/Spectral.h
functionspace/Spectral.cc
)

list( APPEND atlas_numerics_srcs
numerics/Method.h
numerics/Method.cc
numerics/Nabla.h
numerics/Nabla.cc
numerics/fvm/Method.h
numerics/fvm/Method.cc
numerics/fvm/Nabla.h
numerics/fvm/Nabla.cc
)
if( ATLAS_HAVE_TRANS )
list( APPEND atlas_numerics_srcs
  trans/Trans.h
  trans/Trans.cc
)
endif()

list( APPEND atlas_interpolation_srcs
interpolation/Vector2D.h
interpolation/Vector3D.h
interpolation/element/Quad3D.cc
interpolation/element/Quad3D.h
interpolation/element/Quad2D.cc
interpolation/element/Quad2D.h
interpolation/element/Triag3D.cc
interpolation/element/Triag3D.h
interpolation/element/Triag2D.cc
interpolation/element/Triag2D.h
interpolation/method/FiniteElement.cc
interpolation/method/FiniteElement.h
interpolation/method/Intersect.cc
interpolation/method/Intersect.h
interpolation/method/KNearestNeighbours.cc
interpolation/method/KNearestNeighbours.h
interpolation/method/KNearestNeighboursBase.cc
interpolation/method/KNearestNeighboursBase.h
interpolation/method/Method.cc
interpolation/method/Method.h
interpolation/method/NearestNeighbour.cc
interpolation/method/NearestNeighbour.h
interpolation/method/PointIndex3.cc
interpolation/method/PointIndex3.h
interpolation/method/PointSet.cc
interpolation/method/PointSet.h
interpolation/method/Ray.cc
interpolation/method/Ray.h
)

list( APPEND atlas_array_srcs
array/Array.h
array/Array.cc
array/ArrayUtil.h
array/ArrayUtil.cc
array/ArrayView_iterator.h
array/ArrayView_impl.h
array/IndexView.h
array/IndexView.cc
array/ArrayView.h
array/ArrayView.cc
array/DataType.h
)

list( APPEND atlas_util_srcs
util/Config.h
util/Config.cc
util/Constants.h
util/Metadata.h
util/Metadata.cc
output/detail/GmshIO.h
output/detail/GmshIO.cc
output/detail/PointCloud.h
output/detail/PointCloud.cc
parallel/omp/omp.h
parallel/omp/omp.cc
parallel/mpi/Buffer.h
parallel/mpi/mpi.h
parallel/mpi/mpi.cc
parallel/Checksum.h
parallel/Checksum.cc
parallel/HaloExchange.h
parallel/HaloExchange.cc
parallel/GatherScatter.h
parallel/GatherScatter.cc
)

list( APPEND atlas_internals_srcs
mesh/detail/AccumulateFacets.h
mesh/detail/AccumulateFacets.cc
util/Bitflags.h
util/Checksum.h
util/Checksum.cc
runtime/Debug.h
util/MicroDeg.h
mesh/IsGhostNode.h
util/LonLatMicroDeg.h
util/CoordinateEnums.h
mesh/detail/PeriodicTransform.h
util/Unique.h
parallel/detail/MPLArrayView.h
)



### atlas c++ library

ecbuild_debug_var(CGAL_LIBRARIES)
ecbuild_debug_var(CGAL_INCLUDE_DIRS)

if( NOT ATLAS_HAVE_TRANS )
  unset( TRANSI_INCLUDE_DIRS )
  unset( TRANSI_LIBRARIES )
endif()

if( ECKIT_HAVE_MPI AND ATLAS_HAVE_MPI AND  EC_OS_NAME MATCHES "macosx" )
  list( APPEND ATLAS_DEFINITIONS BUG_ECKIT_166 )
  link_libraries(${MPI_CXX_LIBRARIES})
endif()

ecbuild_add_library( TARGET atlas

  AUTO_VERSION

  LINKER_LANGUAGE CXX

  INSTALL_HEADERS ALL

  HEADER_DESTINATION  ${INSTALL_INCLUDE_DIR}/atlas

  SOURCES
    ${atlas_srcs}
    ${atlas_array_srcs}
    ${atlas_grid_srcs}
    ${atlas_mesh_srcs}
    ${atlas_field_srcs}
    ${atlas_functionspace_srcs}
    ${atlas_interpolation_srcs}
    ${atlas_numerics_srcs}
    ${atlas_output_srcs}
    ${atlas_util_srcs}
    ${atlas_internals_srcs}
    ${CMAKE_CURRENT_BINARY_DIR}/library/git_sha1.h
    ${CMAKE_CURRENT_BINARY_DIR}/library/defines.h

  PRIVATE_INCLUDES
    "${CGAL_INCLUDE_DIRS}"
    "${TRANSI_INCLUDE_DIRS}"
    "${MPI_CXX_INCLUDE_DIRS}"

  LIBS
    eckit_geometry
    eckit_linalg
    eckit_mpi
    eckit_option
    "${CGAL_LIBRARIES}"
    "${TRANSI_LIBRARIES}"

  DEFINITIONS
    ${ATLAS_DEFINITIONS}
 )
