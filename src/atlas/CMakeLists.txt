# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

### config headers

ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( atlas_config.h.in   atlas_config.h  )
configure_file( atlas_version.h.in  atlas_version.h )
configure_file( atlas_defines.h.in  atlas_defines.h   )

install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/atlas_version.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( atlas_version.cc.in atlas_version.cc )

### sources


list( APPEND atlas_srcs
atlas.h
atlas.cc
)

list( APPEND atlas_grid_srcs
grid/BoundBox.h
grid/BoundBox.cc
grid/Domain.h
grid/Domain.cc
grid/GaussianGrid.h
grid/GaussianGrid.cc
grid/GaussianLatitudes.h
grid/GaussianLatitudes.cc
grid/Grid.h
grid/Grid.cc
grid/GridDistribution.h
grid/GridDistribution.cc
grid/grids.h
grid/grids.cc
grid/LocalGrid.h
grid/LocalGrid.cc
grid/LonLatGrid.h
grid/LonLatGrid.cc
grid/OctahedralReducedGaussianGrid.h
grid/OctahedralReducedGaussianGrid.cc
grid/ReducedGaussianGrid.h
grid/ReducedGaussianGrid.cc
grid/ReducedGrid.h
grid/ReducedGrid.cc
grid/ReducedLonLatGrid.h
grid/ReducedLonLatGrid.cc
grid/RotatedGrid.h
grid/RotatedGrid.cc
grid/Unstructured.h
grid/Unstructured.cc
grid/predefined/gausslat/gausslat.h
grid/predefined/gausslat/gausslat.cc
grid/predefined/gausslat/N16.cc    # TL31
grid/predefined/gausslat/N24.cc    # TL47
grid/predefined/gausslat/N32.cc    # TL63
grid/predefined/gausslat/N48.cc    # TL95
grid/predefined/gausslat/N64.cc    # TL127
grid/predefined/gausslat/N80.cc    # TL159
grid/predefined/gausslat/N96.cc    # TL191
grid/predefined/gausslat/N128.cc   # TL255
grid/predefined/gausslat/N160.cc   # TL319
grid/predefined/gausslat/N200.cc   # TL399
grid/predefined/gausslat/N256.cc   # TL511
grid/predefined/gausslat/N320.cc   # TL639
grid/predefined/gausslat/N400.cc   # TL799
grid/predefined/gausslat/N512.cc   # TL1023
grid/predefined/gausslat/N576.cc   # TL1151
grid/predefined/gausslat/N640.cc   # TL1279
grid/predefined/gausslat/N800.cc   # TL1599
grid/predefined/gausslat/N1024.cc  # TL2047
grid/predefined/gausslat/N1280.cc  # TL2559
grid/predefined/gausslat/N1600.cc  # TL3199
grid/predefined/gausslat/N2000.cc  # TL3999
grid/predefined/gausslat/N4000.cc  # TL7999
grid/predefined/gausslat/N8000.cc  # TL15999
grid/predefined/rgg/rgg.h
grid/predefined/rgg/N16.cc    # TL31
grid/predefined/rgg/N24.cc    # TL47
grid/predefined/rgg/N32.cc    # TL63
grid/predefined/rgg/N48.cc    # TL95
grid/predefined/rgg/N64.cc    # TL127
grid/predefined/rgg/N80.cc    # TL159
grid/predefined/rgg/N96.cc    # TL191
grid/predefined/rgg/N128.cc   # TL255
grid/predefined/rgg/N160.cc   # TL319
grid/predefined/rgg/N200.cc   # TL399
grid/predefined/rgg/N256.cc   # TL511
grid/predefined/rgg/N320.cc   # TL639
grid/predefined/rgg/N400.cc   # TL799
grid/predefined/rgg/N512.cc   # TL1023
grid/predefined/rgg/N576.cc   # TL1151
grid/predefined/rgg/N640.cc   # TL1279
grid/predefined/rgg/N800.cc   # TL1599
grid/predefined/rgg/N1024.cc  # TL2047
grid/predefined/rgg/N1280.cc  # TL2559
grid/predefined/rgg/N1600.cc  # TL3199
grid/predefined/rgg/N2000.cc  # TL3999
grid/predefined/rgg/N4000.cc  # TL7999
grid/predefined/rgg/N8000.cc  # TL15999
grid/partitioners/EqualRegionsPartitioner.h
grid/partitioners/EqualRegionsPartitioner.cc
grid/partitioners/Partitioner.h
grid/partitioners/Partitioner.cc
)
if( ATLAS_HAVE_TRANS )
list( APPEND atlas_grid_srcs
  grid/partitioners/TransPartitioner.h
  grid/partitioners/TransPartitioner.cc
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
mesh/actions/AddVirtualNodes.h
mesh/actions/AddVirtualNodes.cc
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
mesh/generators/Delaunay.h
mesh/generators/Delaunay.cc
mesh/generators/MeshGenerator.cc
mesh/generators/MeshGenerator.h
mesh/generators/ReducedGridMeshGenerator.h
mesh/generators/ReducedGridMeshGenerator.cc
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
)

list( APPEND atlas_functionspace_srcs
functionspace/Edges.h
functionspace/Edges.cc
functionspace/FunctionSpace.h
functionspace/FunctionSpace.cc
functionspace/Nodes.h
functionspace/Nodes.cc
functionspace/NodesInterface.h
functionspace/NodesInterface.cc
functionspace/ReducedGridPoint.h
functionspace/ReducedGridPoint.cc
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
  numerics/trans/Trans.h
  numerics/trans/Trans.cc
)
endif()

list( APPEND atlas_interpolation_srcs
interpolation/Intersect.h
interpolation/Intersect.cc
interpolation/Quad3D.h
interpolation/Quad3D.cc
interpolation/Ray.h
interpolation/Ray.cc
interpolation/Triag3D.h
interpolation/Triag3D.cc
interpolation/Vector3D.h
interpolation/PointIndex3.h
interpolation/PointIndex3.cc
interpolation/PointSet.h
interpolation/PointSet.cc
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
)

list( APPEND atlas_util_srcs
util/Config.h
util/Config.cc
util/Constants.h
util/DataType.h
util/Metadata.h
util/Metadata.cc
util/io/Gmsh.h
util/io/Gmsh.cc
util/io/PointCloud.h
util/io/PointCloud.cc
util/parallel/atlas_omp.h
util/parallel/atlas_omp.cc
util/parallel/mpi/mpi.h
util/parallel/mpi/mpi.cc
util/parallel/mpi/Collectives.h
util/parallel/mpl/Checksum.h
util/parallel/mpl/Checksum.cc
util/parallel/mpl/HaloExchange.h
util/parallel/mpl/HaloExchange.cc
util/parallel/mpl/GatherScatter.h
util/parallel/mpl/GatherScatter.cc
util/parallel/mpl/MPLArrayView.h
util/runtime/Behavior.h
util/runtime/Behavior.cc
util/runtime/LogFormat.h
util/runtime/LogFormat.cc
util/runtime/ErrorHandling.cc
util/runtime/ErrorHandling.cc
)

list( APPEND atlas_internals_srcs
internals/AccumulateFaces.h
internals/AccumulateFaces.cc
internals/Bitflags.h
internals/Checksum.h
internals/Checksum.cc
internals/Debug.h
internals/Functions.h
internals/IsGhost.h
internals/LonLatMicroDeg.h
internals/ObjectRegistry.h
internals/Parameters.h
internals/PeriodicTransform.h
internals/Unique.h
)



### atlas c++ library

ecbuild_debug_var(CGAL_LIBRARIES)
ecbuild_debug_var(CGAL_INCLUDE_DIRS)

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
    ${atlas_util_srcs}
    ${atlas_internals_srcs}
    ${CMAKE_CURRENT_BINARY_DIR}/atlas_version.h
    ${CMAKE_CURRENT_BINARY_DIR}/atlas_version.cc
    ${CMAKE_CURRENT_BINARY_DIR}/atlas_config.h
    ${CMAKE_CURRENT_BINARY_DIR}/atlas_defines.h

  PRIVATE_INCLUDES
    "${CGAL_INCLUDE_DIRS}"
    "${TRANSI_INCLUDE_DIRS}"

  LIBS
    eckit_geometry
    "${ECKIT_MPI_LIBRARIES}"
    "${CGAL_LIBRARIES}"
    "${TRANSI_LIBRARIES}" )
