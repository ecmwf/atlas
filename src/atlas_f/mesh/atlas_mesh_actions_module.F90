! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_mesh_actions_module

implicit none

public :: atlas_build_parallel_fields
public :: atlas_build_nodes_parallel_fields
public :: atlas_build_edges_parallel_fields
public :: atlas_build_periodic_boundaries
public :: atlas_build_halo
public :: atlas_build_edges
public :: atlas_build_pole_edges
public :: atlas_build_node_to_cell_connectivity
public :: atlas_build_node_to_edge_connectivity
public :: atlas_build_median_dual_mesh
public :: atlas_write_load_balance_report

! =============================================================================
CONTAINS
! =============================================================================

subroutine atlas_build_parallel_fields(mesh)
  use atlas_BuildParallelFields_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_parallel_fields(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_parallel_fields

subroutine atlas_build_nodes_parallel_fields(nodes)
  use atlas_BuildParallelFields_c_binding
  use atlas_mesh_Nodes_module, only: atlas_mesh_Nodes
  type(atlas_mesh_Nodes), intent(inout) :: nodes
  call atlas__build_nodes_parallel_fields(nodes%CPTR_PGIBUG_A)
end subroutine atlas_build_nodes_parallel_fields

subroutine atlas_renumber_nodes_glb_idx(nodes)
  use atlas_BuildParallelFields_c_binding
  use atlas_mesh_Nodes_module, only: atlas_mesh_Nodes
  type(atlas_mesh_Nodes), intent(inout) :: nodes
  call atlas__renumber_nodes_glb_idx(nodes%CPTR_PGIBUG_A)
end subroutine atlas_renumber_nodes_glb_idx

subroutine atlas_build_edges_parallel_fields(mesh)
  use atlas_BuildParallelFields_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_edges_parallel_fields(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_edges_parallel_fields

subroutine atlas_build_periodic_boundaries(mesh)
  use atlas_BuildPeriodicBoundaries_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_periodic_boundaries(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_periodic_boundaries

subroutine atlas_build_halo(mesh,nelems)
  use atlas_BuildHalo_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in) :: nelems
  call atlas__build_halo(mesh%CPTR_PGIBUG_A,nelems)
end subroutine atlas_build_halo

subroutine atlas_build_edges(mesh)
  use atlas_BuildEdges_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_edges(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_edges

subroutine atlas_build_pole_edges(mesh)
  use atlas_BuildEdges_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_pole_edges(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_pole_edges

subroutine atlas_build_node_to_cell_connectivity(mesh)
  use atlas_BuildNode2CellConnectivity_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_node_to_cell_connectivity(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_node_to_cell_connectivity

subroutine atlas_build_node_to_edge_connectivity(mesh)
  use atlas_BuildEdges_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_node_to_edge_connectivity(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_node_to_edge_connectivity

subroutine atlas_build_median_dual_mesh(mesh)
  use atlas_BuildDualMesh_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_median_dual_mesh(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_median_dual_mesh

subroutine atlas_build_centroid_dual_mesh(mesh)
  use atlas_BuildDualMesh_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_centroid_dual_mesh(mesh%CPTR_PGIBUG_A)
end subroutine atlas_build_centroid_dual_mesh

subroutine atlas_write_load_balance_report(mesh,filename)
  use fckit_c_interop_module, only : c_str
  use atlas_WriteLoadBalanceReport_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  type(atlas_Mesh), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_load_balance_report(mesh%CPTR_PGIBUG_A,c_str(filename))
end subroutine atlas_write_load_balance_report

! -----------------------------------------------------------------------------

end module atlas_mesh_actions_module
