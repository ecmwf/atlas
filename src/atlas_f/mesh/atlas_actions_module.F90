! (C) Copyright 2013-2015 ECMWF.

#include "atlas_f/atlas_f_defines.h"

module atlas_actions_module

use atlas_c_interop, only : c_str
use atlas_FunctionSpace_module, only: &
    & atlas_FunctionSpace
use atlas_Field_module, only: &
    & atlas_Field
use atlas_FieldSet_module, only: &
    & atlas_FieldSet
use atlas_Nodes_module, only: &
    & atlas_Nodes
use atlas_Mesh_module, only: &
    & atlas_Mesh
use atlas_Grid_module, only: &
    & atlas_Grid
use atlas_GridDistribution_module, only: &
    & atlas_GridDistribution
use atlas_functionspace_NodeColumns_module, only: &
    & atlas_functionspace_NodeColumns

implicit none

public :: atlas_generate_mesh
public :: atlas_read_gmsh
public :: atlas_write_gmsh
public :: atlas_write_gmsh_field
public :: atlas_write_gmsh_fieldset
public :: atlas_build_parallel_fields
public :: atlas_build_nodes_parallel_fields
public :: atlas_build_edges_parallel_fields
public :: atlas_build_periodic_boundaries
public :: atlas_build_halo
public :: atlas_build_edges
public :: atlas_build_pole_edges
public :: atlas_build_node_to_edge_connectivity
public :: atlas_build_median_dual_mesh
public :: atlas_write_load_balance_report

public

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

!private

!------------------------------------------------------------------------------

interface atlas_generate_mesh
  module procedure atlas_generate_mesh
  module procedure atlas_generate_mesh_with_distribution
end interface atlas_generate_mesh

! =============================================================================
CONTAINS
! =============================================================================

function atlas_read_gmsh(filename) result(mesh)
  use atlas_gmsh_c_binding
  character(len=*), intent(in) :: filename
  type(atlas_Mesh) :: mesh
  mesh = atlas_Mesh( atlas__read_gmsh(c_str(filename)) )
  call mesh%return()
end function atlas_read_gmsh

subroutine atlas_write_gmsh(mesh,filename)
  use atlas_gmsh_c_binding
  type(atlas_Mesh), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_gmsh_mesh(mesh%c_ptr(),c_str(filename))
end subroutine atlas_write_gmsh

subroutine atlas_write_gmsh_field(field,function_space,filename,mode)
  use atlas_gmsh_c_binding
  type(atlas_Field), intent(in) :: field
  type(atlas_functionspace_NodeColumns), intent(in) :: function_space
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_field(field%c_ptr(),function_space%c_ptr(),c_str(filename),mode)
  else
    call atlas__write_gmsh_field(field%c_ptr(),function_space%c_ptr(),c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_field

subroutine atlas_write_gmsh_fieldset(fieldset,function_space,filename,mode)
  use atlas_gmsh_c_binding
  type(atlas_FieldSet), intent(in) :: fieldset
  type(atlas_functionspace_NodeColumns), intent(in) :: function_space
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_fieldset(fieldset%c_ptr(),function_space%c_ptr(),c_str(filename),mode)
  else
    call atlas__write_gmsh_fieldset(fieldset%c_ptr(),function_space%c_ptr(),c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_fieldset

subroutine atlas_build_parallel_fields(mesh)
  use atlas_BuildParallelFields_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_parallel_fields(mesh%c_ptr())
end subroutine atlas_build_parallel_fields

subroutine atlas_build_nodes_parallel_fields(nodes)
  use atlas_BuildParallelFields_c_binding
  type(atlas_Nodes), intent(inout) :: nodes
  call atlas__build_nodes_parallel_fields(nodes%c_ptr())
end subroutine atlas_build_nodes_parallel_fields

subroutine atlas_renumber_nodes_glb_idx(nodes)
  use atlas_BuildParallelFields_c_binding
  type(atlas_Nodes), intent(inout) :: nodes
  call atlas__renumber_nodes_glb_idx(nodes%c_ptr())
end subroutine atlas_renumber_nodes_glb_idx

subroutine atlas_build_edges_parallel_fields(mesh)
  use atlas_BuildParallelFields_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_edges_parallel_fields(mesh%c_ptr())
end subroutine atlas_build_edges_parallel_fields

subroutine atlas_build_periodic_boundaries(mesh)
  use atlas_BuildPeriodicBoundaries_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_periodic_boundaries(mesh%c_ptr())
end subroutine atlas_build_periodic_boundaries

subroutine atlas_build_halo(mesh,nelems)
  use atlas_BuildHalo_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in) :: nelems
  call atlas__build_halo(mesh%c_ptr(),nelems)
end subroutine atlas_build_halo

subroutine atlas_build_edges(mesh)
  use atlas_BuildEdges_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_edges(mesh%c_ptr())
end subroutine atlas_build_edges

subroutine atlas_build_pole_edges(mesh)
  use atlas_BuildEdges_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_pole_edges(mesh%c_ptr())
end subroutine atlas_build_pole_edges

subroutine atlas_build_node_to_edge_connectivity(mesh)
  use atlas_BuildEdges_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_node_to_edge_connectivity(mesh%c_ptr())
end subroutine atlas_build_node_to_edge_connectivity

subroutine atlas_build_median_dual_mesh(mesh)
  use atlas_BuildDualMesh_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_median_dual_mesh(mesh%c_ptr())
end subroutine atlas_build_median_dual_mesh

#if !DEPRECATE_OLD_FUNCTIONSPACE
subroutine atlas_build_centroid_dual_mesh(mesh)
  use atlas_BuildDualMesh_c_binding
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_centroid_dual_mesh(mesh%c_ptr())
end subroutine atlas_build_centroid_dual_mesh
#endif

subroutine atlas_write_load_balance_report(mesh,filename)
  use atlas_WriteLoadBalanceReport_c_binding
  type(atlas_Mesh), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_load_balance_report(mesh%c_ptr(),c_str(filename))
end subroutine atlas_write_load_balance_report

function atlas_generate_mesh(grid) result(mesh)
  use atlas_GenerateMesh_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_Grid) :: grid
  mesh = atlas_Mesh( atlas__generate_mesh(grid%c_ptr()) )
  call mesh%return()
end function atlas_generate_mesh

function atlas_generate_mesh_with_distribution(grid,distribution) result(mesh)
  use atlas_GenerateMesh_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_Grid) :: grid
  type(atlas_GridDistribution) :: distribution
  mesh = atlas_Mesh( atlas__generate_mesh_with_distribution(grid%c_ptr(),distribution%c_ptr()) )
  call mesh%return()
end function atlas_generate_mesh_with_distribution

! -----------------------------------------------------------------------------

end module atlas_actions_module
