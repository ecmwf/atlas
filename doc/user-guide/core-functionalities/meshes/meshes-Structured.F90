program main
use atlas_module
implicit none
type(atlas_Grid)          :: grid
type(atlas_Mesh)          :: mesh
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Output)        :: gmsh_2d, gmsh_3d

call atlas_init()

! Generate mesh
meshgenerator = atlas_meshgenerator_Structured()
grid = atlas_StructuredGrid( "O32" )
mesh = meshgenerator%generate(grid)

gmsh_2d = atlas_output_Gmsh("mesh2d.msh")
gmsh_3d = atlas_output_Gmsh("mesh3d.msh",coordinates="xyz")

! Write mesh
call gmsh_2d%write(mesh)
call gmsh_3d%write(mesh)

! Cleanup
call grid%final()
call mesh%final()
call gmsh_2d%final()
call gmsh_3d%final()
call meshgenerator%final()

call atlas_finalize()

end program main
