program main

use atlas_module

character(len=1024)       :: gridID, visualize
type(atlas_ReducedGrid)   :: reducedGrid
type(atlas_mesh)          :: mesh
type(atlas_meshgenerator) :: meshgenerator
call atlas_init()

call atlas_resource("--grid", "N32", gridID)

reducedGrid = atlas_ReducedGrid(gridID)

meshgenerator = atlas_reducedgridmeshgenerator()
mesh = meshgenerator%generate(reducedGrid)
call atlas_write_gmsh(mesh, "mesh.msh")

call reducedGrid%final()
call meshgenerator%final()
call mesh%final()

call atlas_finalize()

end program main

