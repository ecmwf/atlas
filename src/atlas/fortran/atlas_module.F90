! (C) Copyright 2013-2014 ECMWF.

#include "atlas/atlas_defines_fortran.h"

module atlas_module

! Purpose :
! -------
!    *atlas* : Low-level API to manipulate meshes in a
!        object-oriented fashion
!
!        Objects of classes defined in this module don't contain
!        any data. They only contain a pointer to a object created
!        in a C++ equivalent object, which takes care of all the
!        memory management. Object member functions give access to the
!        real data.

! Classes :
! -------
!   Mesh_type
!   FunctionSpace_type
!   Field_type
!   FieldSet_type
!   Metadata_type
!   HaloExchange_type

! Interfaces :
! ----------

! Author :
! ------
!   20-Nov-2013 Willem Deconinck    *ECMWF*

!------------------------------------------------------------------------------
use, intrinsic :: iso_c_binding
use atlas_mpl_module
use atlas_C_interop
use atlas_atlas_c_binding
use atlas_field_c_binding
use atlas_fieldset_c_binding
use atlas_functionspace_c_binding
use atlas_mesh_c_binding
use atlas_metadata_c_binding
use atlas_haloexchange_c_binding
use atlas_gatherscatter_c_binding
use atlas_rgg_c_binding
use atlas_checksum_c_binding
use atlas_gmsh_c_binding
use atlas_BuildPeriodicBoundaries_c_binding
use atlas_BuildEdges_c_binding
use atlas_BuildDualMesh_c_binding
use atlas_BuildParallelFields_c_binding
use atlas_BuildHalo_c_binding
use atlas_GenerateMesh_c_binding
use atlas_WriteLoadBalanceReport_c_binding
use atlas_atlas_logging_c_binding
implicit none

! ----------------------------------------------------
! ENUM FieldType
integer, private, parameter :: KIND_INT32  = -4
integer, private, parameter :: KIND_REAL32 =  4
integer, private, parameter :: KIND_REAL64 =  8
! ----------------------------------------------------

integer, private, parameter :: FIELD_NB_VARS = -1
integer, private, parameter :: wp = c_double ! working precision

#include "atlas_module_Logging_i.f"
#include "atlas_module_HaloExchange_i.f"
#include "atlas_module_GatherScatter_i.f"
#include "atlas_module_Grid_i.f"
#include "atlas_module_Checksum_i.f"
#include "atlas_module_Mesh_i.f"
#include "atlas_module_FunctionSpace_i.f"
#include "atlas_module_Field_i.f"
#include "atlas_module_FieldSet_i.f"
#include "atlas_module_Metadata_i.f"

INTERFACE delete

! Purpose :
! -------
!   *delete* : Common interface to properly call the destructor
!              of class objects

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

! -----------------------------------------------------------------------------
  module procedure Mesh__delete
  module procedure FunctionSpace__delete
  module procedure FieldSet__delete
  module procedure HaloExchange__delete
  module procedure Metadata__delete
end interface delete

!------------------------------------------------------------------------------

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

!------------------------------------------------------------------------------

! Logger singleton
TYPE(Logger_t) :: logger

! =============================================================================
CONTAINS
! =============================================================================

subroutine atlas_init()
  integer, save :: argc
  type(c_ptr), save :: argv(15)
  call get_c_arguments(argc,argv)
  call atlas__atlas_init(argc,argv)
  !call atlas__atlas_init_noargs()
end subroutine

subroutine atlas_finalize()
  call atlas__atlas_finalize()
end subroutine


integer function real_kind(kind)
  integer :: kind
  if (kind == c_double) then
    real_kind = KIND_REAL64
  else if (kind == c_float) then
    real_kind = KIND_REAL32
  else
    write(0,*) "Unsupported kind"
    write(0,*) 'call abort()'
  end if
end function

integer function integer_kind(kind)
  integer, optional :: kind
  integer_kind = KIND_INT32
  if ( present(kind) ) then
    if (kind == c_int) then
      integer_kind = KIND_INT32
    else
      write(0,*) "Unsupported kind"
      write(0,*) 'call abort()'
    end if
  end if
end function

#include "atlas_module_Logging_c.f"
#include "atlas_module_HaloExchange_c.f"
#include "atlas_module_GatherScatter_c.f"
#include "atlas_module_Grid_c.f"
#include "atlas_module_Checksum_c.f"
#include "atlas_module_Mesh_c.f"
#include "atlas_module_FunctionSpace_c.f"
#include "atlas_module_Field_c.f"
#include "atlas_module_FieldSet_c.f"
#include "atlas_module_Metadata_c.f"



! -----------------------------------------------------------------------------

function eckit_version()
  character(len=5) :: eckit_version
  eckit_version = c_to_f_string_cptr(atlas__eckit_version())
end function eckit_version

function eckit_git_sha1()
  character(len=40) :: eckit_git_sha1
  eckit_git_sha1 = c_to_f_string_cptr(atlas__eckit_git_sha1())
end function eckit_git_sha1

function atlas_version()
  character(len=5) :: atlas_version
  atlas_version = c_to_f_string_cptr(atlas__atlas_version())
end function atlas_version

function atlas_git_sha1()
  character(len=40) :: atlas_git_sha1
  atlas_git_sha1 = c_to_f_string_cptr(atlas__atlas_git_sha1())
end function atlas_git_sha1

function atlas_read_gmsh(filename) result(mesh)
  character(len=*), intent(in) :: filename
  type(Mesh_type) :: mesh
  mesh%private%object = atlas__read_gmsh(c_str(filename))
end function atlas_read_gmsh

subroutine atlas_write_gmsh(mesh,filename)
  type(Mesh_type), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_gmsh_mesh(mesh%private%object,c_str(filename))
end subroutine atlas_write_gmsh

subroutine atlas_write_gmsh_field(field,filename,mode)
  type(Field_type), intent(in) :: field
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_field(field%private%object,c_str(filename),mode)
  else
    call atlas__write_gmsh_field(field%private%object,c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_field

subroutine atlas_write_gmsh_fieldset(fieldset,filename,mode)
  type(FieldSet_type), intent(in) :: fieldset
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_fieldset(fieldset%private%object,c_str(filename),mode)
  else
    call atlas__write_gmsh_fieldset(fieldset%private%object,c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_fieldset

subroutine atlas_build_parallel_fields(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_parallel_fields(mesh%private%object)
end subroutine atlas_build_parallel_fields

subroutine atlas_build_nodes_parallel_fields(nodes)
  type(FunctionSpace_type), intent(inout) :: nodes
  call atlas__build_nodes_parallel_fields(nodes%private%object)
end subroutine atlas_build_nodes_parallel_fields

subroutine atlas_renumber_nodes_glb_idx(nodes)
  type(FunctionSpace_type), intent(inout) :: nodes
  call atlas__renumber_nodes_glb_idx(nodes%private%object)
end subroutine atlas_renumber_nodes_glb_idx

subroutine atlas_build_edges_parallel_fields(edges, nodes)
  type(FunctionSpace_type), intent(inout) :: edges, nodes
  call atlas__build_edges_parallel_fields(edges%private%object,nodes%private%object)
end subroutine atlas_build_edges_parallel_fields

subroutine atlas_build_periodic_boundaries(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_periodic_boundaries(mesh%private%object)
end subroutine atlas_build_periodic_boundaries

subroutine atlas_build_halo(mesh,nelems)
  type(Mesh_type), intent(inout) :: mesh
  integer, intent(in) :: nelems
  call atlas__build_halo(mesh%private%object,nelems)
end subroutine atlas_build_halo

subroutine atlas_build_edges(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_edges(mesh%private%object)
end subroutine atlas_build_edges

subroutine atlas_build_pole_edges(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_pole_edges(mesh%private%object)
end subroutine atlas_build_pole_edges

subroutine atlas_build_node_to_edge_connectivity(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_node_to_edge_connectivity(mesh%private%object)
end subroutine atlas_build_node_to_edge_connectivity

subroutine atlas_build_median_dual_mesh(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_median_dual_mesh(mesh%private%object)
end subroutine atlas_build_median_dual_mesh

subroutine atlas_build_centroid_dual_mesh(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_centroid_dual_mesh(mesh%private%object)
end subroutine atlas_build_centroid_dual_mesh

subroutine atlas_write_load_balance_report(mesh,filename)
  type(Mesh_type), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_load_balance_report(mesh%private%object,c_str(filename))
end subroutine atlas_write_load_balance_report

subroutine atlas_generate_reduced_gaussian_grid(mesh,identifier)
  type(Mesh_type), intent(inout) :: mesh
  character(len=*), intent(in) :: identifier
  mesh%private%object = atlas__generate_reduced_gaussian_grid(c_str(identifier))
end subroutine atlas_generate_reduced_gaussian_grid

subroutine atlas_generate_full_gaussian_grid(mesh,nlon,nlat)
  type(Mesh_type), intent(inout) :: mesh
  integer, intent(in) :: nlon, nlat
  mesh%private%object = atlas__generate_full_gaussian_grid(nlon,nlat)
end subroutine atlas_generate_full_gaussian_grid

subroutine atlas_generate_latlon_grid(mesh,nlon,nlat)
  type(Mesh_type), intent(inout) :: mesh
  integer, intent(in) :: nlon, nlat
  mesh%private%object = atlas__generate_latlon_grid(nlon,nlat)
end subroutine atlas_generate_latlon_grid

subroutine atlas_generate_custom_reduced_gaussian_grid(mesh,nlon)
  type(Mesh_type), intent(inout) :: mesh
  integer, intent(in) :: nlon(:)
  mesh%private%object = atlas__generate_custom_reduced_gaussian_grid(nlon,size(nlon))
end subroutine atlas_generate_custom_reduced_gaussian_grid

function atlas_generate_mesh(rgg) result(mesh)
  type(Mesh_type) :: mesh
  type(ReducedGG_type) :: rgg
  mesh%private%object = atlas__generate_mesh(rgg%private%object)
end function atlas_generate_mesh


! -----------------------------------------------------------------------------

end module atlas_module
