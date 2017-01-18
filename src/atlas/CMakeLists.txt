# (C) Copyright 1996-2016 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

### config headers

ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( internals/atlas_config.h.in   internals/atlas_config.h  )
configure_file( internals/atlas_version.h.in  internals/atlas_version.h )
configure_file( internals/atlas_defines.h.in  internals/atlas_defines.h   )

install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_defines.h
  ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_version.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas/internals )

configure_file( internals/atlas_version.cc.in internals/atlas_version.cc )

### sources

list( APPEND atlas_srcs
atlas.h
atlas.cc
)

list( APPEND atlas_grid_srcs
grid/CustomStructured.cc
grid/CustomStructured.h
grid/Domain.cc
grid/Domain.h
grid/Grid.cc
grid/Grid.h
grid/GridDistribution.cc
grid/GridDistribution.h
grid/Structured.cc
grid/Structured.h
grid/Unstructured.cc
grid/Unstructured.h
grid/gaussian/ClassicGaussian.cc
grid/gaussian/ClassicGaussian.h
grid/gaussian/Gaussian.cc
grid/gaussian/Gaussian.h
grid/gaussian/OctahedralGaussian.cc
grid/gaussian/OctahedralGaussian.h
grid/gaussian/ReducedGaussian.cc
grid/gaussian/ReducedGaussian.h
grid/gaussian/RegularGaussian.cc
grid/gaussian/RegularGaussian.h
grid/gaussian/classic/N.cc
grid/gaussian/classic/N.h
grid/gaussian/classic/N1024.cc  # TL2047
grid/gaussian/classic/N128.cc   # TL255
grid/gaussian/classic/N1280.cc  # TL2559
grid/gaussian/classic/N16.cc    # TL31
grid/gaussian/classic/N160.cc   # TL319
grid/gaussian/classic/N1600.cc  # TL3199
grid/gaussian/classic/N200.cc   # TL399
grid/gaussian/classic/N2000.cc  # TL3999
grid/gaussian/classic/N24.cc    # TL47
grid/gaussian/classic/N256.cc   # TL511
grid/gaussian/classic/N32.cc    # TL63
grid/gaussian/classic/N320.cc   # TL639
grid/gaussian/classic/N400.cc   # TL799
grid/gaussian/classic/N4000.cc  # TL7999
grid/gaussian/classic/N48.cc    # TL95
grid/gaussian/classic/N512.cc   # TL1023
grid/gaussian/classic/N576.cc   # TL1151
grid/gaussian/classic/N64.cc    # TL127
grid/gaussian/classic/N640.cc   # TL1279
grid/gaussian/classic/N80.cc    # TL159
grid/gaussian/classic/N800.cc   # TL1599
grid/gaussian/classic/N8000.cc  # TL15999
grid/gaussian/classic/N96.cc    # TL191
grid/gaussian/classic/PointsPerLatitude.cc
grid/gaussian/classic/PointsPerLatitude.h
grid/gaussian/latitudes/Latitudes.cc
grid/gaussian/latitudes/Latitudes.h
grid/gaussian/latitudes/N.cc
grid/gaussian/latitudes/N.h
grid/gaussian/latitudes/N1024.cc  # TL2047
grid/gaussian/latitudes/N128.cc   # TL255
grid/gaussian/latitudes/N1280.cc  # TL2559
grid/gaussian/latitudes/N16.cc    # TL31
grid/gaussian/latitudes/N160.cc   # TL319
grid/gaussian/latitudes/N1600.cc  # TL3199
grid/gaussian/latitudes/N200.cc   # TL399
grid/gaussian/latitudes/N2000.cc  # TL3999
grid/gaussian/latitudes/N24.cc    # TL47
grid/gaussian/latitudes/N256.cc   # TL511
grid/gaussian/latitudes/N32.cc    # TL63
grid/gaussian/latitudes/N320.cc   # TL639
grid/gaussian/latitudes/N400.cc   # TL799
grid/gaussian/latitudes/N4000.cc  # TL7999
grid/gaussian/latitudes/N48.cc    # TL95
grid/gaussian/latitudes/N512.cc   # TL1023
grid/gaussian/latitudes/N576.cc   # TL1151
grid/gaussian/latitudes/N64.cc    # TL127
grid/gaussian/latitudes/N640.cc   # TL1279
grid/gaussian/latitudes/N80.cc    # TL159
grid/gaussian/latitudes/N800.cc   # TL1599
grid/gaussian/latitudes/N8000.cc  # TL15999
grid/gaussian/latitudes/N96.cc    # TL191
grid/grids.cc
grid/grids.h
grid/lonlat/LonLat.cc
grid/lonlat/LonLat.h
grid/lonlat/ReducedLonLat.cc
grid/lonlat/ReducedLonLat.h
grid/lonlat/RegularLonLat.cc
grid/lonlat/RegularLonLat.h
grid/lonlat/ShiftedLat.cc
grid/lonlat/ShiftedLat.h
grid/lonlat/ShiftedLon.cc
grid/lonlat/ShiftedLon.h
grid/lonlat/ShiftedLonLat.cc
grid/lonlat/ShiftedLonLat.h
grid/partitioners/EqualRegionsPartitioner.cc
grid/partitioners/EqualRegionsPartitioner.h
grid/partitioners/Partitioner.cc
grid/partitioners/Partitioner.h
)
if( ATLAS_HAVE_TRANS )
list( APPEND atlas_grid_srcs
  grid/partitioners/TransPartitioner.h
  grid/partitioners/TransPartitioner.cc
)
endif()

list( APPEND atlas_mesh_srcs
mesh/Connectivity.cc
mesh/Connectivity.h
mesh/ElementType.cc
mesh/ElementType.h
mesh/Elements.cc
mesh/Elements.h
mesh/Halo.cc
mesh/Halo.h
mesh/HybridElements.cc
mesh/HybridElements.h
mesh/Mesh.cc
mesh/Mesh.h
mesh/Nodes.cc
mesh/Nodes.h
mesh/actions/AddVirtualNodes.cc
mesh/actions/AddVirtualNodes.h
mesh/actions/BuildCellCentres.cc
mesh/actions/BuildCellCentres.h
mesh/actions/BuildConvexHull3D.cc
mesh/actions/BuildConvexHull3D.h
mesh/actions/BuildDualMesh.cc
mesh/actions/BuildDualMesh.h
mesh/actions/BuildEdges.cc
mesh/actions/BuildEdges.h
mesh/actions/BuildHalo.cc
mesh/actions/BuildHalo.h
mesh/actions/BuildParallelFields.cc
mesh/actions/BuildParallelFields.h
mesh/actions/BuildPeriodicBoundaries.cc
mesh/actions/BuildPeriodicBoundaries.h
mesh/actions/BuildStatistics.cc
mesh/actions/BuildStatistics.h
mesh/actions/BuildXYZField.cc
mesh/actions/BuildXYZField.h
mesh/actions/GenerateMesh.cc
mesh/actions/GenerateMesh.h
mesh/actions/WriteLoadBalanceReport.cc
mesh/actions/WriteLoadBalanceReport.h
mesh/generators/Delaunay.cc
mesh/generators/Delaunay.h
mesh/generators/MeshGenerator.cc
mesh/generators/MeshGenerator.h
mesh/generators/Structured.cc
mesh/generators/Structured.h

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
array/array.h
array/array_fwd.h
array/Array.h
array/ArrayIdx.h
array/ArrayLayout.h
array/ArrayShape.h
array/ArraySpec.cc
array/ArraySpec.h
array/ArrayStrides.h
array/ArrayUtil.h
array/ArrayView.h
array/DataType.h
array/IndexView.h
array/LocalView.cc
array/LocalView.h
array/StorageView.h
array/Vector.h
)
if( ATLAS_HAVE_GRIDTOOLS_STORAGE )
list( APPEND atlas_array_srcs
array/gridtools/GridToolsArray.cc
array/gridtools/GridToolsArrayHelpers.h
array/gridtools/GridToolsArrayView.cc
array/gridtools/GridToolsArrayView.h
array/gridtools/GridToolsDataStore.h
array/gridtools/GridToolsIndexView.cc
array/gridtools/GridToolsIndexView.h
array/gridtools/GridToolsMakeView.cc
array/gridtools/GridToolsMakeView.h
array/gridtools/GridToolsStorageView.h
array/gridtools/GridToolsTraits.h
)
else()
list( APPEND atlas_array_srcs
array/native/NativeArray.cc
array/native/NativeArrayView.cc
array/native/NativeArrayView.h
array/native/NativeDataStore.h
array/native/NativeIndexView.cc
array/native/NativeIndexView.h
array/native/NativeMakeView.cc
array/native/NativeStorageView.h
)
endif()

list( APPEND atlas_util_srcs
parallel/Checksum.cc
parallel/Checksum.h
parallel/GatherScatter.cc
parallel/GatherScatter.h
parallel/HaloExchange.cc
parallel/HaloExchange.h
parallel/mpi/Collectives.h
parallel/mpi/mpi.cc
parallel/mpi/mpi.h
parallel/omp/omp.cc
parallel/omp/omp.h
runtime/Behavior.cc
runtime/Behavior.h
runtime/ErrorHandling.cc
runtime/ErrorHandling.cc
runtime/LogFormat.cc
runtime/LogFormat.h
util/Config.cc
util/Config.h
util/Constants.h
util/Metadata.cc
util/Metadata.h
util/io/Gmsh.cc
util/io/Gmsh.h
util/io/PointCloud.cc
util/io/PointCloud.h

)

list( APPEND atlas_internals_srcs
internals/AccumulateFaces.cc
internals/AccumulateFaces.h
internals/AtlasTool.cc
internals/AtlasTool.h
internals/Bitflags.h
internals/Checksum.cc
internals/Checksum.h
internals/Functions.h
internals/IsGhost.h
internals/LonLatMicroDeg.h
internals/Parameters.h
internals/PeriodicTransform.h
internals/Unique.h
)



### atlas c++ library

ecbuild_debug_var(CGAL_LIBRARIES)
ecbuild_debug_var(CGAL_INCLUDE_DIRS)

if( NOT ATLAS_HAVE_TRANS )
  unset( TRANSI_INCLUDE_DIRS )
  unset( TRANSI_LIBRARIES )
endif()

if( SHORTCUT_COMPILATION )
  unset( atlas_srcs )
  unset( atlas_grid_srcs )
  unset( atlas_mesh_srcs )
  unset( atlas_field_srcs )
  unset( atlas_functionspace_srcs )
  unset( atlas_interpolation_srcs )
  unset( atlas_numerics_srcs )
  unset( atlas_output_srcs )
  unset( atlas_util_srcs )
  unset( atlas_internals_srcs )
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
    ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_version.h
    ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_version.cc
    ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_config.h
    ${CMAKE_CURRENT_BINARY_DIR}/internals/atlas_defines.h

  PRIVATE_INCLUDES
    "${CGAL_INCLUDE_DIRS}"
    "${TRANSI_INCLUDE_DIRS}"
  LIBS
    eckit_geometry eckit_option
    "${ECKIT_MPI_LIBRARIES}"
    "${CGAL_LIBRARIES}"
    "${TRANSI_LIBRARIES}"
    "${CUDA_LIBRARIES}" )
