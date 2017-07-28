program atlas_interpolation_fortran
use atlas_module
implicit none

type(atlas_Grid) :: source_grid
type(atlas_Grid) :: target_grid

type(atlas_Mesh) :: source_mesh
type(atlas_Mesh) :: target_mesh

type(atlas_functionspace_NodeColumns) :: source_fs
type(atlas_functionspace_NodeColumns) :: target_fs

type(atlas_Field) :: source_field
type(atlas_Field) :: target_field

type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Partitioner) :: target_partitioner
type(atlas_GridDistribution) :: target_distribution

type(atlas_Output) :: gmsh

call atlas_library%initialise()

meshgenerator = atlas_MeshGenerator()

! Generate source mesh
source_grid = atlas_Grid("O20")
source_mesh = meshgenerator%generate(source_grid)

! Generate target mesh based on domain decomposition of source mesh
target_grid = atlas_Grid("L80")
target_partitioner  = atlas_MatchingMeshPartitioner(source_mesh)
target_distribution = target_partitioner%partition(target_grid)
target_mesh = meshgenerator%generate(target_grid,target_distribution)

! Output target mesh
gmsh = atlas_output_Gmsh("target.msh")
call gmsh%write(target_mesh)



source_fs = atlas_functionspace_NodeColumns(source_mesh,halo=1)
target_fs = atlas_functionspace_NodeColumns(target_mesh,halo=0)

source_field = source_fs%create_field("source",atlas_real(8))
target_field = target_fs%create_field("target",atlas_real(8))



! cleanup
call gmsh%final()
call target_partitioner%final()
call target_distribution%final()
call source_field%final()
call target_field%final()
call meshgenerator%final()
call source_fs%final()
call target_fs%final()
call source_mesh%final()
call target_mesh%final()
call source_grid%final()
call target_grid%final()

call atlas_library%finalise()
contains
end program
