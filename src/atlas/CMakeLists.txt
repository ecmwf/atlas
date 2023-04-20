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
library.h
library/config.h
library/FloatingPointExceptions.h
library/FloatingPointExceptions.cc
library/Library.h
library/Library.cc
library/Plugin.h
library/Plugin.cc
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
parallel/omp/copy.h
parallel/omp/fill.h
parallel/omp/sort.h
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
projection/Jacobian.h
projection/detail/CubedSphereEquiAnglProjection.cc
projection/detail/CubedSphereEquiAnglProjection.h
projection/detail/CubedSphereEquiDistProjection.cc
projection/detail/CubedSphereEquiDistProjection.h
projection/detail/CubedSphereProjectionBase.cc
projection/detail/CubedSphereProjectionBase.h
projection/detail/LambertAzimuthalEqualAreaProjection.cc
projection/detail/LambertAzimuthalEqualAreaProjection.h
projection/detail/LambertConformalConicProjection.cc
projection/detail/LambertConformalConicProjection.h
projection/detail/LonLatProjection.cc
projection/detail/LonLatProjection.h
projection/detail/MercatorProjection.cc
projection/detail/MercatorProjection.h
projection/detail/ProjectionFactory.cc
projection/detail/ProjectionFactory.h
projection/detail/ProjectionImpl.cc
projection/detail/ProjectionImpl.h
projection/detail/ProjectionUtilities.h
projection/detail/SchmidtProjection.cc
projection/detail/SchmidtProjection.h
projection/detail/VariableResolutionProjection.cc
projection/detail/VariableResolutionProjection.h

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
grid/CubedSphereGrid.cc
grid/CubedSphereGrid.h
grid/Grid.cc
grid/Grid.h
grid/SpecRegistry.cc
grid/SpecRegistry.h
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
grid/StructuredPartitionPolygon.h
grid/StructuredPartitionPolygon.cc
grid/Tiles.h
grid/Tiles.cc

grid/detail/grid/CubedSphere.cc
grid/detail/grid/CubedSphere.h
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
grid/detail/grid/RegionalVariableResolution.h
grid/detail/grid/RegionalVariableResolution.cc
grid/detail/grid/Healpix.h
grid/detail/grid/Healpix.cc

grid/detail/distribution/DistributionImpl.h
grid/detail/distribution/DistributionImpl.cc
grid/detail/distribution/DistributionArray.cc
grid/detail/distribution/DistributionArray.h
grid/detail/distribution/DistributionFunction.cc
grid/detail/distribution/DistributionFunction.h

grid/detail/distribution/BandsDistribution.cc
grid/detail/distribution/BandsDistribution.h
grid/detail/distribution/SerialDistribution.cc
grid/detail/distribution/SerialDistribution.h

grid/detail/tiles/Tiles.cc
grid/detail/tiles/Tiles.h
grid/detail/tiles/TilesFactory.h
grid/detail/tiles/TilesFactory.cc
grid/detail/tiles/FV3Tiles.cc
grid/detail/tiles/FV3Tiles.h
grid/detail/tiles/LFRicTiles.cc
grid/detail/tiles/LFRicTiles.h

grid/detail/vertical/VerticalInterface.h
grid/detail/vertical/VerticalInterface.cc

grid/detail/partitioner/BandsPartitioner.cc
grid/detail/partitioner/BandsPartitioner.h
grid/detail/partitioner/CheckerboardPartitioner.cc
grid/detail/partitioner/CheckerboardPartitioner.h
grid/detail/partitioner/CubedSpherePartitioner.cc
grid/detail/partitioner/CubedSpherePartitioner.h
grid/detail/partitioner/EqualBandsPartitioner.cc
grid/detail/partitioner/EqualBandsPartitioner.h
grid/detail/partitioner/EqualRegionsPartitioner.cc
grid/detail/partitioner/EqualRegionsPartitioner.h
grid/detail/partitioner/MatchingMeshPartitioner.h
grid/detail/partitioner/MatchingMeshPartitioner.cc
grid/detail/partitioner/MatchingMeshPartitionerBruteForce.cc
grid/detail/partitioner/MatchingMeshPartitionerBruteForce.h
grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.cc
grid/detail/partitioner/MatchingMeshPartitionerCubedSphere.h
grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.cc
grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.h
grid/detail/partitioner/MatchingMeshPartitionerSphericalPolygon.cc
grid/detail/partitioner/MatchingMeshPartitionerSphericalPolygon.h
grid/detail/partitioner/MatchingFunctionSpacePartitioner.h
grid/detail/partitioner/MatchingFunctionSpacePartitioner.cc
grid/detail/partitioner/MatchingFunctionSpacePartitionerLonLatPolygon.cc
grid/detail/partitioner/MatchingFunctionSpacePartitionerLonLatPolygon.h
grid/detail/partitioner/Partitioner.cc
grid/detail/partitioner/Partitioner.h
grid/detail/partitioner/RegularBandsPartitioner.cc
grid/detail/partitioner/RegularBandsPartitioner.h
grid/detail/partitioner/SerialPartitioner.cc
grid/detail/partitioner/SerialPartitioner.h

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
if( atlas_HAVE_TRANS )
list( APPEND atlas_grid_srcs
  grid/detail/partitioner/TransPartitioner.h
  grid/detail/partitioner/TransPartitioner.cc
)
endif()
if( atlas_HAVE_PROJ )
list( APPEND atlas_grid_srcs
  projection/detail/ProjProjection.cc
  projection/detail/ProjProjection.h
)
endif()


# Append CGAL_COMPLE_FLAGS only to this file ( see ATLAS-193 )
if( CGAL_COMPILE_FLAGS )
  set_source_files_properties( mesh/actions/BuildConvexHull3D.cc
                               PROPERTIES COMPILE_FLAGS "${CGAL_COMPILE_FLAGS}" )
endif()

list( APPEND atlas_mesh_srcs
mesh.h
mesh/BuildMeshFromConnectivities.cc
mesh/BuildMeshFromConnectivities.h
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
mesh/actions/Build2DCellCentres.cc
mesh/actions/Build2DCellCentres.h
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
mesh/actions/Reorder.h
mesh/actions/Reorder.cc
mesh/actions/ReorderHilbert.h
mesh/actions/ReorderHilbert.cc
mesh/actions/ReorderReverseCuthillMckee.h
mesh/actions/ReorderReverseCuthillMckee.cc

meshgenerator.h
meshgenerator/MeshGenerator.cc
meshgenerator/MeshGenerator.h
meshgenerator/detail/CubedSphereMeshGenerator.h
meshgenerator/detail/CubedSphereMeshGenerator.cc
meshgenerator/detail/CubedSphereDualMeshGenerator.h
meshgenerator/detail/CubedSphereDualMeshGenerator.cc
meshgenerator/detail/NodalCubedSphereMeshGenerator.h
meshgenerator/detail/NodalCubedSphereMeshGenerator.cc
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
meshgenerator/detail/HealpixMeshGenerator.h
meshgenerator/detail/HealpixMeshGenerator.cc
meshgenerator/detail/cubedsphere/CubedSphereUtility.h
meshgenerator/detail/cubedsphere/CubedSphereUtility.cc
)

list( APPEND atlas_output_srcs
output/Output.h
output/Output.cc
output/Gmsh.h
output/Gmsh.cc
output/detail/GmshIO.cc
output/detail/GmshIO.h
output/detail/GmshImpl.cc
output/detail/GmshImpl.h
output/detail/GmshInterface.cc
output/detail/GmshInterface.h
output/detail/PointCloudIO.cc
output/detail/PointCloudIO.h

)

list( APPEND atlas_field_srcs
field.h
field/Field.cc
field/Field.h
field/FieldCreator.cc
field/FieldCreator.h
field/FieldCreatorArraySpec.cc
field/FieldCreatorArraySpec.h
field/FieldCreatorIFS.cc
field/FieldCreatorIFS.h
field/FieldSet.cc
field/FieldSet.h
field/MissingValue.cc
field/MissingValue.h
field/State.cc
field/State.h
field/detail/FieldImpl.cc
field/detail/FieldImpl.h
field/detail/FieldInterface.cc
field/detail/FieldInterface.h
field/detail/MissingValue.cc
field/detail/MissingValue.h
)

list( APPEND atlas_functionspace_srcs
functionspace.h
functionspace/BlockStructuredColumns.h
functionspace/BlockStructuredColumns.cc
functionspace/CellColumns.h
functionspace/CellColumns.cc
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
functionspace/CubedSphereColumns.h
functionspace/CubedSphereColumns.cc
functionspace/detail/BlockStructuredColumns.h
functionspace/detail/BlockStructuredColumns.cc
functionspace/detail/BlockStructuredColumnsInterface.h
functionspace/detail/BlockStructuredColumnsInterface.cc
functionspace/detail/FunctionSpaceImpl.h
functionspace/detail/FunctionSpaceImpl.cc
functionspace/detail/FunctionSpaceInterface.h
functionspace/detail/FunctionSpaceInterface.cc
functionspace/detail/NodeColumnsInterface.h
functionspace/detail/NodeColumnsInterface.cc
functionspace/detail/NodeColumns_FieldStatistics.cc
functionspace/detail/SpectralInterface.h
functionspace/detail/SpectralInterface.cc
functionspace/detail/StructuredColumns.h
functionspace/detail/StructuredColumns.cc
functionspace/detail/StructuredColumnsInterface.h
functionspace/detail/StructuredColumnsInterface.cc
functionspace/detail/StructuredColumns_setup.cc
functionspace/detail/StructuredColumns_create_remote_index.cc
functionspace/detail/PointCloudInterface.h
functionspace/detail/PointCloudInterface.cc
functionspace/detail/CubedSphereStructure.h
functionspace/detail/CubedSphereStructure.cc

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
if( atlas_HAVE_TRANS )
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
interpolation/Cache.cc
interpolation/Cache.h
interpolation/Interpolation.cc
interpolation/Interpolation.h
interpolation/NonLinear.cc
interpolation/NonLinear.h
interpolation/Vector2D.cc
interpolation/Vector2D.h
interpolation/Vector3D.cc
interpolation/Vector3D.h
interpolation/element/Quad2D.h
interpolation/element/Quad2D.cc
interpolation/element/Quad3D.cc
interpolation/element/Quad3D.h
interpolation/element/Triag2D.cc
interpolation/element/Triag2D.h
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
interpolation/method/PointIndex2.cc
interpolation/method/PointIndex2.h
interpolation/method/PointSet.cc
interpolation/method/PointSet.h
interpolation/method/Ray.cc
interpolation/method/Ray.h
interpolation/method/cubedsphere/CellFinder.cc
interpolation/method/cubedsphere/CellFinder.h
interpolation/method/cubedsphere/CubedSphereBilinear.cc
interpolation/method/cubedsphere/CubedSphereBilinear.h
interpolation/method/knn/GridBox.cc
interpolation/method/knn/GridBox.h
interpolation/method/knn/GridBoxAverage.cc
interpolation/method/knn/GridBoxAverage.h
interpolation/method/knn/GridBoxMaximum.cc
interpolation/method/knn/GridBoxMaximum.h
interpolation/method/knn/GridBoxMethod.cc
interpolation/method/knn/GridBoxMethod.h
interpolation/method/knn/KNearestNeighbours.cc
interpolation/method/knn/KNearestNeighbours.h
interpolation/method/knn/KNearestNeighboursBase.cc
interpolation/method/knn/KNearestNeighboursBase.h
interpolation/method/knn/NearestNeighbour.cc
interpolation/method/knn/NearestNeighbour.h
interpolation/method/structured/Cubic2D.cc
interpolation/method/structured/Cubic2D.h
interpolation/method/structured/Cubic3D.cc
interpolation/method/structured/Cubic3D.h
interpolation/method/structured/Linear2D.cc
interpolation/method/structured/Linear2D.h
interpolation/method/structured/Linear3D.cc
interpolation/method/structured/Linear3D.h
interpolation/method/structured/QuasiCubic2D.cc
interpolation/method/structured/QuasiCubic2D.h
interpolation/method/structured/QuasiCubic3D.cc
interpolation/method/structured/QuasiCubic3D.h
interpolation/method/structured/StructuredInterpolation2D.h
interpolation/method/structured/StructuredInterpolation2D.tcc
interpolation/method/structured/StructuredInterpolation3D.h
interpolation/method/structured/StructuredInterpolation3D.tcc
interpolation/method/structured/kernels/Cubic3DKernel.h
interpolation/method/structured/kernels/CubicHorizontalKernel.h
interpolation/method/structured/kernels/CubicVerticalKernel.h
interpolation/method/structured/kernels/Linear3DKernel.h
interpolation/method/structured/kernels/LinearHorizontalKernel.h
interpolation/method/structured/kernels/LinearVerticalKernel.h
interpolation/method/structured/kernels/QuasiCubic3DKernel.cc
interpolation/method/structured/kernels/QuasiCubic3DKernel.h
interpolation/method/structured/kernels/QuasiCubicHorizontalKernel.h
interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.cc
interpolation/method/unstructured/ConservativeSphericalPolygonInterpolation.h
interpolation/method/unstructured/FiniteElement.cc
interpolation/method/unstructured/FiniteElement.h
interpolation/method/unstructured/UnstructuredBilinearLonLat.cc
interpolation/method/unstructured/UnstructuredBilinearLonLat.h
interpolation/nonlinear/Missing.cc
interpolation/nonlinear/Missing.h
interpolation/nonlinear/NonLinear.cc
interpolation/nonlinear/NonLinear.h

linalg/Indexing.h
linalg/Introspection.h
linalg/View.h
linalg/sparse.h
linalg/sparse/Backend.h
linalg/sparse/Backend.cc
linalg/sparse/SparseMatrixMultiply.h
linalg/sparse/SparseMatrixMultiply.tcc
linalg/sparse/SparseMatrixMultiply_EckitLinalg.h
linalg/sparse/SparseMatrixMultiply_EckitLinalg.cc
linalg/sparse/SparseMatrixMultiply_OpenMP.h
linalg/sparse/SparseMatrixMultiply_OpenMP.cc
linalg/dense.h
linalg/dense/Backend.h
linalg/dense/Backend.cc
linalg/dense/MatrixMultiply.h
linalg/dense/MatrixMultiply.tcc
linalg/dense/MatrixMultiply_EckitLinalg.h
linalg/dense/MatrixMultiply_EckitLinalg.cc
)


list (APPEND atlas_redistribution_srcs
redistribution/Redistribution.h
redistribution/Redistribution.cc
redistribution/detail/RedistributionImpl.h
redistribution/detail/RedistributionImpl.cc
redistribution/detail/RedistributionImplFactory.h
redistribution/detail/RedistributionImplFactory.cc
redistribution/detail/RedistributeGeneric.h
redistribution/detail/RedistributeGeneric.cc
redistribution/detail/RedistributeStructuredColumns.h
redistribution/detail/RedistributeStructuredColumns.cc
)

list( APPEND atlas_array_srcs
array.h
array_fwd.h
array/Array.h
array/ArrayDataStore.cc
array/ArrayDataStore.h
array/ArrayIdx.h
array/ArrayLayout.h
array/ArrayShape.h
array/ArraySpec.cc
array/ArraySpec.h
array/ArrayStrides.h
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
array/helpers/ArrayCopier.h
array/helpers/ArrayForEach.h
#array/Table.h
#array/Table.cc
#array/TableView.h
#array/TableView.cc
)
if( atlas_HAVE_GRIDTOOLS_STORAGE )
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
parallel/Checksum.cc
parallel/Checksum.h
parallel/GatherScatter.cc
parallel/GatherScatter.h
parallel/HaloExchange.cc
parallel/HaloExchange.h
parallel/HaloAdjointExchangeImpl.h
parallel/HaloExchangeImpl.h
parallel/mpi/Buffer.h
runtime/Exception.cc
runtime/Exception.h
util/Config.cc
util/Config.h
util/Constants.h
util/ConvexSphericalPolygon.cc
util/ConvexSphericalPolygon.h
util/Earth.h
util/GaussianLatitudes.cc
util/GaussianLatitudes.h
util/Geometry.cc
util/Geometry.h
util/KDTree.cc
util/KDTree.h
util/PolygonXY.cc
util/PolygonXY.h
util/Metadata.cc
util/Metadata.h
util/Point.cc
util/Point.h
util/Polygon.cc
util/Polygon.h
util/Rotation.cc
util/Rotation.h
util/Registry.h
util/SphericalPolygon.cc
util/SphericalPolygon.h
util/UnitSphere.h
util/vector.h
util/VectorOfAbstract.h
util/detail/BlackMagic.h
util/detail/Cache.h
util/detail/Debug.h
util/detail/KDTree.h
util/function/MDPI_functions.h
util/function/MDPI_functions.cc
util/function/SolidBodyRotation.h
util/function/SolidBodyRotation.cc
util/function/SphericalHarmonic.h
util/function/SphericalHarmonic.cc
util/function/VortexRollup.h
util/function/VortexRollup.cc
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
util/Topology.h
util/Unique.h
util/Unique.cc
util/Allocate.h
util/Allocate.cc
)

list( APPEND atlas_io_srcs
    io/atlas-io.h
    io/Data.cc
    io/Data.h
    io/detail/Base64.cc
    io/detail/Base64.h
    io/detail/Checksum.h
    io/detail/Checksum.cc
    io/detail/DataInfo.h
    io/detail/DataType.cc
    io/detail/DataType.h
    io/detail/Decoder.cc
    io/detail/Decoder.h
    io/detail/Defaults.h
    io/detail/Encoder.cc
    io/detail/Encoder.h
    io/detail/Endian.h
    io/detail/Link.cc
    io/detail/Link.h
    io/detail/ParsedRecord.h
    io/detail/RecordInfo.h
    io/detail/RecordSections.h
    io/detail/Reference.h
    io/detail/sfinae.h
    io/detail/StaticAssert.h
    io/detail/tag.h
    io/detail/Time.cc
    io/detail/Time.h
    io/detail/Type.h
    io/detail/TypeTraits.h
    io/detail/Version.h
    io/Exceptions.cc
    io/Exceptions.h
    io/FileStream.cc
    io/FileStream.h
    io/Metadata.cc
    io/Metadata.h
    io/print/TableFormat.cc
    io/print/TableFormat.h
    io/print/JSONFormat.cc
    io/print/JSONFormat.h
    io/print/Bytes.cc
    io/print/Bytes.h
    io/ReadRequest.cc
    io/ReadRequest.h
    io/Record.cc
    io/Record.h
    io/RecordItem.cc
    io/RecordItem.h
    io/RecordItemReader.cc
    io/RecordItemReader.h
    io/RecordPrinter.cc
    io/RecordPrinter.h
    io/RecordReader.cc
    io/RecordReader.h
    io/RecordWriter.cc
    io/RecordWriter.h
    io/Session.cc
    io/Session.h
    io/Stream.cc
    io/Stream.h
    io/types/array.h
    io/types/array/ArrayMetadata.cc
    io/types/array/ArrayMetadata.h
    io/types/array/ArrayReference.cc
    io/types/array/ArrayReference.h
    io/types/array/adaptors/StdVectorAdaptor.h
    io/types/string.h
    io/types/scalar.h
    io/types/scalar.cc
)


list( APPEND atlas_io_adaptor_srcs
  io/ArrayAdaptor.cc
  io/ArrayAdaptor.h
  io/VectorAdaptor.h
)


### atlas c++ library

if( NOT atlas_HAVE_TRANS )
  unset( TRANSI_INCLUDE_DIRS )
  unset( TRANSI_LIBRARIES )
endif()

if( SHORTCUT_COMPILATION )
  unset( atlas_grid_srcs )
  unset( atlas_mesh_srcs )
  unset( atlas_field_srcs )
  unset( atlas_functionspace_srcs )
  unset( atlas_interpolation_srcs )
  unset( atlas_redistribution_srcs )
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
  ${atlas_redistribution_srcs}
  ${atlas_numerics_srcs}
  ${atlas_output_srcs}
  ${atlas_util_srcs}
  ${atlas_io_adaptor_srcs}
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

#ecbuild_add_library( TARGET atlas_io

#    INSTALL_HEADERS    ALL
#    HEADER_DESTINATION include/atlas_io
#    SOURCES            ${atlas_io_srcs}
#    PUBLIC_LIBS        eckit
#    PUBLIC_INCLUDES
#      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
#      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>


#)


ecbuild_add_library( TARGET atlas

  AUTO_VERSION

  INSTALL_HEADERS ALL

  HEADER_DESTINATION  include/atlas

  SOURCES
    ${source_list}

  PRIVATE_LIBS
    $<${atlas_HAVE_FORTRAN}:fckit>
    $<${atlas_HAVE_TRANS}:transi>
    $<${atlas_HAVE_ACC}:atlas_acc_support>
    ${CGAL_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${PROJ_LIBRARIES}

  PUBLIC_LIBS
    eckit
    eckit_geometry
    eckit_linalg
    eckit_maths
    eckit_mpi
    eckit_option
    atlas_io
    $<${atlas_HAVE_EIGEN}:Eigen3::Eigen>
    $<${atlas_HAVE_OMP_CXX}:OpenMP::OpenMP_CXX>
    $<${atlas_HAVE_GRIDTOOLS_STORAGE}:GridTools::gridtools>

  PRIVATE_INCLUDES
    ${CGAL_INCLUDE_DIRS}
    ${FFTW_INCLUDES}
    ${PROJ_INCLUDE_DIRS}

  PUBLIC_INCLUDES
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
    $<INSTALL_INTERFACE:include>

)

target_compile_features( atlas PUBLIC cxx_std_17 )
