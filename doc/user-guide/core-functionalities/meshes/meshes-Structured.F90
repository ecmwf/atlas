program main

use atlas_module

character(len=1024)       :: gridID, visualize
type(atlas_grid_Structured)   :: Structured
type(atlas_Mesh)          :: mesh
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Output) :: gmsh
call atlas_init()

call atlas_resource("--grid", "N32", gridID)

Structured = atlas_grid_Structured(gridID)

meshgenerator = atlas_meshgenerator_Structured()
mesh = meshgenerator%generate(Structured)
gmsh = atlas_output_Gmsh("mesh.msh")
call gmsh%write(mesh)

call Structured%final()
call meshgenerator%final()
call mesh%final()
call gmsh%final()
call atlas_finalize()

end program main
