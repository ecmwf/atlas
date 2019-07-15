# (C) Copyright 2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

### config headers

ecbuild_generate_config_headers( DESTINATION ${INSTALL_INCLUDE_DIR}/atlas )

configure_file( library/defines.h.in   library/defines.h   )
configure_file( library/version.h.in   library/version.h   )
configure_file( library/git_sha1.h.in  library/git_sha1.h  )

install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/library/defines.h
  ${CMAKE_CURRENT_BINARY_DIR}/library/version.h
  ${CMAKE_CURRENT_BINARY_DIR}/library/git_sha1.h
  DESTINATION
  ${INSTALL_INCLUDE_DIR}/atlas/library
)

### sources

list( APPEND atlas_srcs
library/config.h
library/Library.h
library/Library.cc
runtime/AtlasTool.h
runtime/AtlasTool.cc
runtime/Log.h
runtime/Log.cc
runtime/Trace.h
runtime/trace/CallStack.h
runtime/trace/CallStack.cc
runtime/trace/CodeLocation.cc
runtime/trace/CodeLocation.h
runtime/trace/TraceT.h
runtime/trace/Nesting.cc
runtime/trace/Nesting.h
runtime/trace/Barriers.cc
runtime/trace/Barriers.h
runtime/trace/Logging.cc
runtime/trace/Logging.h
runtime/trace/Timings.h
runtime/trace/Timings.cc
parallel/mpi/mpi.cc
parallel/mpi/mpi.h
parallel/omp/omp.cc
parallel/omp/omp.h
)

list( APPEND atlas_grid_srcs

option.h
option/Options.h
option/Options.cc
option/TransOptions.h
option/TransOptions.cc


projection.h
projection/Projection.cc
projection/Projection.h
projection/detail/LambertAzimuthalEqualAreaProjection.cc
projection/detail/LambertAzimuthalEqualAreaProjection.h
projection/detail/LambertProjection.cc
projection/detail/LambertProjection.h
projection/detail/LonLatProjection.cc
projection/detail/LonLatProjection.h
projection/detail/MercatorProjection.cc
projection/detail/MercatorProjection.h
projection/detail/ProjectionFactory.cc
projection/detail/ProjectionFactory.h
projection/detail/ProjectionImpl.cc
projection/detail/ProjectionImpl.h
projection/detail/SchmidtProjection.cc
projection/detail/SchmidtProjection.h

domain.h
domain/Domain.cc
domain/Domain.h
domain/detail/Domain.cc
domain/detail/Domain.h
domain/detail/RectangularDomain.h
domain/detail/RectangularDomain.cc
domain/detail/EmptyDomain.h
domain/detail/EmptyDomain.cc
domain/detail/GlobalDomain.h
domain/detail/GlobalDomain.cc
domain/detail/ZonalBandDomain.h
domain/detail/ZonalBandDomain.cc
domain/detail/DomainFactory.h
domain/detail/DomainFactory.cc

grid.h
grid/Grid.cc
grid/Grid.h
grid/StructuredGrid.cc
grid/StructuredGrid.h
grid/UnstructuredGrid.cc
grid/UnstructuredGrid.h

grid/Distribution.cc
grid/Distribution.h
grid/Spacing.cc
grid/Spacing.h
grid/Partitioner.h
grid/Partitioner.cc
grid/Iterator.h
grid/Iterator.cc
grid/Vertical.h
grid/Vertical.cc
grid/Stencil.h
grid/StencilComputer.h
grid/StencilComputer.cc

grid/detail/grid/GridBuilder.h
grid/detail/grid/GridBuilder.cc
grid/detail/grid/GridFactory.h
grid/detail/grid/GridFactory.cc
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

grid/detail/distribution/DistributionImpl.h
grid/detail/distribution/DistributionImpl.cc

grid/detail/vertical/VerticalInterface.h
grid/detail/vertical/VerticalInterface.cc

grid/detail/partitioner/CheckerboardPartitioner.cc
grid/detail/partitioner/CheckerboardPartitioner.h
grid/detail/partitioner/EqualRegionsPartitioner.cc
grid/detail/partitioner/EqualRegionsPartitioner.h
grid/detail/partitioner/MatchingMeshPartitioner.h
grid/detail/partitioner/MatchingMeshPartitioner.cc
grid/detail/partitioner/MatchingMeshPartitionerBruteForce.cc
grid/detail/partitioner/MatchingMeshPartitionerBruteForce.h
grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.cc
grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.h
grid/detail/partitioner/MatchingMeshPartitionerSphericalPolygon.cc
grid/detail/partitioner/MatchingMeshPartitionerSphericalPolygon.h
grid/detail/partitioner/Partitioner.cc
grid/detail/partitioner/Partitioner.h

grid/detail/spacing/Spacing.cc
grid/detail/spacing/Spacing.h
grid/detail/spacing/SpacingFactory.cc
grid/detail/spacing/SpacingFactory.h
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

# Append CGAL_COMPLE_FLAGS only to this file ( see ATLAS-193 )
if( CGAL_COMPILE_FLAGS )
  set_source_files_properties( mesh/actions/BuildConvexHull3D.cc
                               PROPERTIES COMPILE_FLAGS "${CGAL_COMPILE_FLAGS}" )
endif()

list( APPEND atlas_mesh_srcs
mesh.h
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
mesh/PartitionPolygon.cc
mesh/PartitionPolygon.h
mesh/detail/MeshImpl.cc
mesh/detail/MeshImpl.h
mesh/detail/MeshIntf.cc
mesh/detail/MeshIntf.h
mesh/detail/PartitionGraph.cc
mesh/detail/PartitionGraph.h

mesh/actions/ExtendNodesGlobal.h
mesh/actions/ExtendNodesGlobal.cc
mesh/actions/BuildDualMesh.h
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
mesh/actions/BuildNode2CellConnectivity.cc
mesh/actions/BuildNode2CellConnectivity.h
mesh/actions/BuildParallelFields.cc
mesh/actions/BuildParallelFields.h
mesh/actions/BuildPeriodicBoundaries.cc
mesh/actions/BuildPeriodicBoundaries.h
mesh/actions/BuildStatistics.cc
mesh/actions/BuildStatistics.h
mesh/actions/BuildXYZField.cc
mesh/actions/BuildXYZField.h
mesh/actions/WriteLoadBalanceReport.cc
mesh/actions/BuildTorusXYZField.h
mesh/actions/BuildTorusXYZField.cc

meshgenerator.h
meshgenerator/MeshGenerator.cc
meshgenerator/MeshGenerator.h
meshgenerator/detail/DelaunayMeshGenerator.h
meshgenerator/detail/DelaunayMeshGenerator.cc
meshgenerator/detail/StructuredMeshGenerator.h
meshgenerator/detail/StructuredMeshGenerator.cc
meshgenerator/detail/RegularMeshGenerator.cc
meshgenerator/detail/RegularMeshGenerator.h
meshgenerator/detail/MeshGeneratorFactory.cc
meshgenerator/detail/MeshGeneratorFactory.h
meshgenerator/detail/MeshGeneratorImpl.cc
meshgenerator/detail/MeshGeneratorImpl.h
meshgenerator/detail/MeshGeneratorInterface.cc
meshgenerator/detail/MeshGeneratorInterface.h
)

list( APPEND atlas_output_srcs
output/Output.h
output/Output.cc
output/Gmsh.h
output/Gmsh.cc
)

list( APPEND atlas_field_srcs
field.h
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
field/detail/FieldImpl.h
field/detail/FieldImpl.cc
field/detail/FieldInterface.h
field/detail/FieldInterface.cc
)

list( APPEND atlas_functionspace_srcs
functionspace.h
functionspace/EdgeColumns.h
functionspace/EdgeColumns.cc
functionspace/FunctionSpace.h
functionspace/FunctionSpace.cc
functionspace/NodeColumns.h
functionspace/NodeColumns.cc
functionspace/StructuredColumns.h
functionspace/StructuredColumns.cc
functionspace/Spectral.h
functionspace/Spectral.cc
functionspace/PointCloud.h
functionspace/PointCloud.cc
functionspace/detail/FunctionSpaceImpl.h
functionspace/detail/FunctionSpaceImpl.cc
functionspace/detail/FunctionSpaceInterface.h
functionspace/detail/FunctionSpaceInterface.cc
functionspace/detail/NodeColumnsInterface.h
functionspace/detail/NodeColumnsInterface.cc
functionspace/detail/NodeColumns_FieldStatistics.cc
functionspace/detail/SpectralInterface.h
functionspace/detail/SpectralInterface.cc
functionspace/detail/StructuredColumnsInterface.h
functionspace/detail/StructuredColumnsInterface.cc
functionspace/detail/StructuredColumns_setup.cc
functionspace/detail/StructuredColumns_create_remote_index.cc
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

trans/Cache.h
trans/Cache.cc
trans/Trans.h
trans/Trans.cc
trans/VorDivToUV.h
trans/VorDivToUV.cc
trans/LegendreCacheCreator.h
trans/LegendreCacheCreator.cc
trans/local/TransLocal.h
trans/local/TransLocal.cc
trans/local/LegendrePolynomials.h
trans/local/LegendrePolynomials.cc
trans/local/VorDivToUVLocal.h
trans/local/VorDivToUVLocal.cc
trans/local/LegendreCacheCreatorLocal.h
trans/local/LegendreCacheCreatorLocal.cc
trans/detail/TransFactory.h
trans/detail/TransFactory.cc
trans/detail/TransImpl.h
trans/detail/TransImpl.cc
trans/detail/TransInterface.h
trans/detail/TransInterface.cc
)
if( ATLAS_HAVE_TRANS )
list( APPEND atlas_numerics_srcs
  trans/ifs/LegendreCacheCreatorIFS.h
  trans/ifs/LegendreCacheCreatorIFS.cc
  trans/ifs/TransIFS.h
  trans/ifs/TransIFS.cc
  trans/ifs/TransIFSNodeColumns.h
  trans/ifs/TransIFSNodeColumns.cc
  trans/ifs/TransIFSStructuredColumns.h
  trans/ifs/TransIFSStructuredColumns.cc
  trans/ifs/VorDivToUVIFS.h
  trans/ifs/VorDivToUVIFS.cc
)
endif()

list( APPEND atlas_interpolation_srcs
interpolation.h
interpolation/Interpolation.h
interpolation/Interpolation.cc
interpolation/Vector2D.cc
interpolation/Vector2D.h
interpolation/Vector3D.cc
interpolation/Vector3D.h
interpolation/element/Quad3D.cc
interpolation/element/Quad3D.h
interpolation/element/Triag3D.cc
interpolation/element/Triag3D.h
interpolation/method/Intersect.cc
interpolation/method/Intersect.h
interpolation/method/Method.cc
interpolation/method/Method.h
interpolation/method/MethodFactory.cc
interpolation/method/MethodFactory.h
interpolation/method/PointIndex3.cc
interpolation/method/PointIndex3.h
interpolation/method/PointSet.cc
interpolation/method/PointSet.h
interpolation/method/Ray.cc
interpolation/method/Ray.h
interpolation/method/fe/FiniteElement.cc
interpolation/method/fe/FiniteElement.h
interpolation/method/knn/KNearestNeighbours.cc
interpolation/method/knn/KNearestNeighbours.h
interpolation/method/knn/KNearestNeighboursBase.cc
interpolation/method/knn/KNearestNeighboursBase.h
interpolation/method/knn/NearestNeighbour.cc
interpolation/method/knn/NearestNeighbour.h
interpolation/method/structured/StructuredInterpolation2D.tcc
interpolation/method/structured/StructuredInterpolation2D.h
interpolation/method/structured/StructuredInterpolation3D.tcc
interpolation/method/structured/StructuredInterpolation3D.h
interpolation/method/structured/kernels/LinearHorizontalKernel.h
interpolation/method/structured/kernels/Linear3DKernel.h
interpolation/method/structured/kernels/LinearVerticalKernel.h
interpolation/method/structured/kernels/CubicHorizontalKernel.h
interpolation/method/structured/kernels/Cubic3DKernel.h
interpolation/method/structured/kernels/CubicVerticalKernel.h
interpolation/method/structured/kernels/QuasiCubicHorizontalKernel.h
interpolation/method/structured/kernels/QuasiCubic3DKernel.cc
interpolation/method/structured/kernels/QuasiCubic3DKernel.h
interpolation/method/structured/Linear2D.cc
interpolation/method/structured/Linear2D.h
interpolation/method/structured/Linear3D.cc
interpolation/method/structured/Linear3D.h
interpolation/method/structured/Cubic2D.cc
interpolation/method/structured/Cubic2D.h
interpolation/method/structured/Cubic3D.cc
interpolation/method/structured/Cubic3D.h
interpolation/method/structured/QuasiCubic2D.cc
interpolation/method/structured/QuasiCubic2D.h
interpolation/method/structured/QuasiCubic3D.cc
interpolation/method/structured/QuasiCubic3D.h
)



list( APPEND atlas_array_srcs
array.h
array_fwd.h
array/Array.h
array/ArrayIdx.h
array/ArrayLayout.h
array/ArrayShape.h
array/ArraySpec.cc
array/ArraySpec.h
array/ArrayStrides.h
array/ArrayUtil.cc
array/ArrayUtil.h
array/ArrayView.h
array/ArrayViewUtil.h
array/ArrayViewDefs.h
array/DataType.cc
array/DataType.h
array/IndexView.h
array/LocalView.cc
array/LocalView.h
array/Range.h
array/Vector.h
array/Vector.cc
array/SVector.h
array/helpers/ArrayInitializer.h
array/helpers/ArrayAssigner.h
array/helpers/ArrayWriter.h
array/helpers/ArraySlicer.h
#array/Table.h
#array/Table.cc
#array/TableView.h
#array/TableView.cc
)
if( ATLAS_HAVE_GRIDTOOLS_STORAGE )
list( APPEND atlas_array_srcs
array/gridtools/GPUClonable.h
array/gridtools/GridToolsArray.cc
array/gridtools/GridToolsArrayHelpers.h
array/gridtools/GridToolsArrayView.cc
array/gridtools/GridToolsArrayView.h
array/gridtools/GridToolsDataStore.h
array/gridtools/GridToolsIndexView.cc
array/gridtools/GridToolsIndexView.h
array/gridtools/GridToolsMakeView.cc
array/gridtools/GridToolsMakeView.h
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
)
endif()

list( APPEND atlas_util_srcs
output/detail/GmshIO.cc
output/detail/GmshIO.h
output/detail/PointCloudIO.cc
output/detail/PointCloudIO.h
parallel/Checksum.cc
parallel/Checksum.h
parallel/GatherScatter.cc
parallel/GatherScatter.h
parallel/HaloExchange.cc
parallel/HaloExchange.h
parallel/HaloExchangeImpl.h
parallel/mpi/Buffer.h
runtime/Exception.cc
runtime/Exception.h
util/Config.cc
util/Config.h
util/Constants.h
util/Earth.h
util/GaussianLatitudes.cc
util/GaussianLatitudes.h
util/LonLatPolygon.cc
util/LonLatPolygon.h
util/Metadata.cc
util/Metadata.h
util/Point.cc
util/Point.h
util/Polygon.cc
util/Polygon.h
util/Rotation.cc
util/Rotation.h
util/SphericalPolygon.cc
util/SphericalPolygon.h
util/UnitSphere.h
util/detail/BlackMagic.h
util/detail/Cache.h
util/detail/Debug.h
)

list( APPEND atlas_internals_srcs
mesh/detail/AccumulateFacets.h
mesh/detail/AccumulateFacets.cc
util/Object.h
util/Object.cc
util/ObjectHandle.h
util/ObjectHandle.cc
util/Factory.h
util/Factory.cc
util/Bitflags.h
util/Checksum.h
util/Checksum.cc
util/MicroDeg.h
mesh/IsGhostNode.h
util/LonLatMicroDeg.h
util/CoordinateEnums.h
util/PeriodicTransform.h
util/Unique.h
util/Unique.cc
util/Allocate.h
util/Allocate.cc
#parallel/detail/MPLArrayView.h
)

### atlas c++ library

if( NOT ATLAS_HAVE_TRANS )
  unset( TRANSI_INCLUDE_DIRS )
  unset( TRANSI_LIBRARIES )
endif()

if( SHORTCUT_COMPILATION )
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

list( APPEND source_list
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
)

if( ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA )
  list( APPEND source_list
    parallel/HaloExchangeCUDA.h
    parallel/HaloExchangeCUDA.cu
  )
endif()

include( atlas_host_device )
atlas_host_device( source_list
  SOURCES
    array/gridtools/GridToolsArrayView.cc
    array/gridtools/GridToolsIndexView.cc
    mesh/Connectivity.cc
)

ecbuild_add_library( TARGET atlas

  AUTO_VERSION

  INSTALL_HEADERS ALL

  HEADER_DESTINATION  include/atlas

  SOURCES
    ${source_list}
)

if( ATLAS_HAVE_FORTRAN )
  target_link_libraries( atlas PRIVATE fckit )
endif()

if( ATLAS_HAVE_TRANS )
  if( TRANSI_VERSION VERSION_LESS 0.6 )
    target_link_libraries( atlas PRIVATE transi )
    target_include_directories( atlas PRIVATE ${TRANSI_INCLUDE_DIRS} )
  else()
    if( MIR_NEEDS_TRANSI )
      target_link_libraries( atlas PUBLIC transi ) # Can become private when MIR doesn't use transi_version()
    else()
      target_link_libraries( atlas PRIVATE transi )
    endif()
  endif()
endif()

if( eckit_VERSION )
    set( eckit_link_public PUBLIC )
else()
    set( eckit_link_public PRIVATE )
endif()

target_link_libraries( atlas ${eckit_link_public}
  eckit
  eckit_geometry
  eckit_linalg
  eckit_mpi
  eckit_option
)
if( ECKIT_INCLUDE_DIRS )
    target_include_directories( atlas PRIVATE ${ECKIT_INCLUDE_DIRS} )
endif()

target_include_directories( atlas PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
    $<INSTALL_INTERFACE:include>
)



if( ATLAS_HAVE_TESSELATION )
  target_link_libraries( atlas PRIVATE ${CGAL_LIBRARIES} )
  target_include_directories( atlas PRIVATE ${CGAL_INCLUDE_DIRS} )
endif()

if( ATLAS_HAVE_FFTW )
  target_link_libraries( atlas PRIVATE ${FFTW_LIBRARIES} )
  target_include_directories( atlas PRIVATE ${FFTW_INCLUDES} )
endif()

if( ATLAS_HAVE_EIGEN AND EIGEN3_INCLUDE_DIR )
  target_include_directories( atlas PUBLIC ${EIGEN3_INCLUDE_DIR} )
endif()

if( ATLAS_HAVE_OMP_CXX )
  target_link_libraries( atlas PRIVATE OpenMP::OpenMP_CXX )
endif()

if( ATLAS_HAVE_GRIDTOOLS_STORAGE )
  target_link_libraries( atlas PUBLIC gridtools::storage )
endif()

if( ATLAS_HAVE_ACC )
  target_link_libraries( atlas PRIVATE atlas_acc_support )
endif()

