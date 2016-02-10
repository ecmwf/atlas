program main

use atlas_module

character(len=1024)             :: gridID, visualize
type(atlas_ReducedGrid)         :: reducedGrid
type(atlas_mesh)                :: mesh
type(atlas_meshgenerator)       :: mg
type(atlas_functionspace_nodes) :: fs_nodes

call atlas_init()

call atlas_resource("--grid", "N32", gridID)
reducedGrid = atlas_ReducedGrid(gridID)

mg = atlas_reducedgridmeshgenerator()
mesh          = mg%generate(reducedGrid)
fs_nodes      = atlas_functionspace_nodes(mesh)

call atlas_write_gmsh(mesh, "mesh.msh")

call reducedGrid%final()
call mesh%final()
call mg%final()
call fs_nodes%final()

call atlas_finalize()

end program main

