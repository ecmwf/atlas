! (C) Copyright 2013-2014 ECMWF.

#include "atlas_f/atlas_f_defines.h"

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
!   atlas_Mesh
!   atlas_FunctionSpace
!   atlas_Field
!   atlas_FieldSet
!   atlas_Metadata
!   atlas_HaloExchange

! Interfaces :
! ----------

! Author :
! ------
!   20-Nov-2013 Willem Deconinck    *ECMWF*

!------------------------------------------------------------------------------
use, intrinsic :: iso_c_binding
use atlas_mpi_module
use atlas_C_interop
use atlas_atlas_c_binding
use atlas_mpi_c_binding
use atlas_Field_c_binding
use atlas_Fieldset_c_binding
use atlas_FunctionSpace_c_binding
use atlas_Mesh_c_binding
use atlas_Metadata_c_binding
use atlas_haloexchange_c_binding
use atlas_gatherscatter_c_binding
use atlas_grids_c_binding
use atlas_reducedgrid_c_binding
use atlas_griddistribution_c_binding
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

private :: object_type
private :: view1d
private :: stride
private :: get_c_arguments
private :: c_to_f_string_str
private :: c_to_f_string_cptr
private :: c_str
private :: c_str_no_trim
private :: resource_get_int32
private :: resource_get_int64
private :: resource_get_real32
private :: resource_get_real64
private :: resource_get_string
private :: resource_set_int32
private :: resource_set_int64
private :: resource_set_real32
private :: resource_set_real64
private :: resource_set_string

! ----------------------------------------------------
! ENUM FieldType
integer, public, parameter :: ATLAS_KIND_INT32  = -4
integer, public, parameter :: ATLAS_KIND_INT64  = -8
integer, public, parameter :: ATLAS_KIND_REAL32 =  4
integer, public, parameter :: ATLAS_KIND_REAL64 =  8
! ----------------------------------------------------

integer, private, parameter :: FIELD_NB_VARS = 2147483647 ! maximum integer value
integer, public, parameter :: ATLAS_FIELD_NB_VARS = FIELD_NB_VARS ! maximum integer value
integer, private, parameter :: wp = c_double ! working precision

#if ATLAS_BITS_GLOBAL == 32
integer, public, parameter :: ATLAS_KIND_GIDX = c_int
#elif ATLAS_BITS_GLOBAL == 64
integer, public, parameter :: ATLAS_KIND_GIDX = c_long
#else
#error ATLAS_BITS_GLOBAL must be either 32 or 64
#endif

#include "atlas_module_Config_i.f"
#include "atlas_module_Logging_i.f"
#include "atlas_module_HaloExchange_i.f"
#include "atlas_module_GatherScatter_i.f"
#include "atlas_module_Grid_i.f"
#include "atlas_module_Checksum_i.f"
#include "atlas_module_Mesh_i.f"
#include "atlas_module_FunctionSpace_i.f"
#include "atlas_module_Field_i.f"
#include "atlas_module_FieldSet_i.f"
#include "atlas_module_JSON_i.f"
#include "atlas_module_Metadata_i.f"
#include "atlas_module_NodesFunctionSpace_i.f"
#include "atlas_module_PathName_i.f"
#include "atlas_module_Error_i.f"
#include "atlas_module_GridDistribution_i.f"
#include "atlas_module_State_i.f"
#include "atlas_module_Trans_i.f"
#include "atlas_module_Value_i.f"

INTERFACE atlas_delete

! Purpose :
! -------
!   *delete* : Common interface to properly call the destructor
!              of class objects

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

! -----------------------------------------------------------------------------
  module procedure atlas_ReducedGrid__delete
  module procedure atlas_Mesh__delete
  module procedure atlas_FieldSet__delete
  module procedure atlas_HaloExchange__delete
  module procedure atlas_Metadata__delete
  module procedure atlas_NodesFunctionSpace__delete
  module procedure atlas_SpectralFunctionSpace__delete
  module procedure atlas_Config__delete
  module procedure atlas_ConfigList__delete
  module procedure atlas_GridDistribution__delete
  module procedure atlas_Trans__delete
  module procedure atlas_Value__delete
  module procedure atlas_Value__array_delete
  module procedure atlas_State__delete
end interface atlas_delete

!------------------------------------------------------------------------------

INTERFACE atlas_resource

! Purpose :
! -------
!   *resource* : Configuration
!
! Author :
! ------
!   20-dec-2014 Willem Deconinck     *ECMWF*
! -----------------------------------------------------------------------------
  module procedure resource_get_int32
  module procedure resource_get_int64
  module procedure resource_get_real32
  module procedure resource_get_real64
  module procedure resource_get_string
end interface atlas_resource

!------------------------------------------------------------------------------

INTERFACE atlas_resource_set

! Purpose :
! -------
!   *resource* : Configuration
!
! Author :
! ------
!   10-june-2015 Willem Deconinck     *ECMWF*
! -----------------------------------------------------------------------------
  module procedure resource_set_int32
  module procedure resource_set_int64
  module procedure resource_set_real32
  module procedure resource_set_real64
  module procedure resource_set_string
end interface atlas_resource_set

!------------------------------------------------------------------------------

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

!------------------------------------------------------------------------------

interface atlas_generate_mesh
  module procedure atlas_generate_mesh
  module procedure atlas_generate_mesh_with_distribution
end interface atlas_generate_mesh


!------------------------------------------------------------------------------

! Logger singleton
TYPE(atlas_Logger) :: atlas_log

! =============================================================================
CONTAINS
! =============================================================================


subroutine resource_get_int32(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_int), intent(in) :: default_value
  integer(c_int), intent(out) :: value
  value = atlas__resource_int( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_int64(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_long), intent(in) :: default_value
  integer(c_long), intent(out) :: value
  value = atlas__resource_long( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_real32(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_float), intent(in) :: default_value
  real(c_float), intent(out) :: value
  value = atlas__resource_float( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_real64(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_double), intent(in) :: default_value
  real(c_double), intent(out) :: value
  value = atlas__resource_double( c_str(resource_str), default_value )
end subroutine

subroutine resource_get_string(resource_str,default_value,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  character(len=*), intent(in) :: default_value
  character(len=*), intent(out) :: value
  type(c_ptr) :: value_c_str
  value_c_str = atlas__resource_string( c_str(resource_str), c_str(default_value) )
  value = c_to_f_string_cptr(value_c_str)
end subroutine

subroutine resource_set_int32(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_int), intent(in) :: value
  call atlas__resource_set_int( c_str(resource_str), value )
end subroutine

subroutine resource_set_int64(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  integer(c_long), intent(in) ::value
  call atlas__resource_set_long( c_str(resource_str), value )
end subroutine

subroutine resource_set_real32(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_float), intent(in) :: value
  call atlas__resource_set_float( c_str(resource_str), value )
end subroutine

subroutine resource_set_real64(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  real(c_double), intent(in) :: value
  call atlas__resource_set_double( c_str(resource_str), value )
end subroutine

subroutine resource_set_string(resource_str,value)
  use atlas_atlas_resource_c_binding
  character(len=*), intent(in) :: resource_str
  character(len=*), intent(in) :: value
  type(c_ptr) :: value_c_str
  call atlas__resource_set_string( c_str(resource_str), c_str(value) )
end subroutine

! =============================================================================


subroutine atlas_init()
  integer, save :: argc
  type(c_ptr), save :: argv(15)
  call get_c_arguments(argc,argv)
  call atlas__atlas_init(argc,argv)
  atlas_log = atlas_Logger()
  !call atlas__atlas_init_noargs()
end subroutine

subroutine atlas_finalize()
  call atlas__atlas_finalize()
end subroutine

integer function atlas_real(kind)
  integer :: kind
  if (kind == c_double) then
    atlas_real = ATLAS_KIND_REAL64
  else if (kind == c_float) then
    atlas_real = ATLAS_KIND_REAL32
  else
    call atlas_abort("Unsupported real kind")
  end if
end function

integer function atlas_integer(kind)
  integer, optional :: kind
  atlas_integer = ATLAS_KIND_INT32
  if ( present(kind) ) then
    if (kind == c_int) then
      atlas_integer = ATLAS_KIND_INT32
    else if (kind == c_long) then
      atlas_integer = ATLAS_KIND_INT64
    else
      call atlas_abort("Unsupported real kind")
    end if
  end if
end function

function atlas_data_type(kind)
  character(len=6) :: atlas_data_type
  integer, intent(in) :: kind
  if( kind == ATLAS_KIND_INT32 ) then
    atlas_data_type = "int32"
  else if( kind == ATLAS_KIND_INT64 ) then
    atlas_data_type = "int64"
  else if( kind == ATLAS_KIND_REAL32 ) then
    atlas_data_type = "real32"
  else if( kind == ATLAS_KIND_REAL64 ) then
    atlas_data_type = "real64"
  else
    call atlas_abort("cannot convert kind to data_type",atlas_code_location(__FILE__,__LINE__))
  endif
end function

#include "atlas_module_Config_c.f"
#include "atlas_module_Logging_c.f"
#include "atlas_module_HaloExchange_c.f"
#include "atlas_module_GatherScatter_c.f"
#include "atlas_module_Grid_c.f"
#include "atlas_module_Checksum_c.f"
#include "atlas_module_Mesh_c.f"
#include "atlas_module_FunctionSpace_c.f"
#include "atlas_module_Field_c.f"
#include "atlas_module_FieldSet_c.f"
#include "atlas_module_JSON_c.f"
#include "atlas_module_Metadata_c.f"
#include "atlas_module_NodesFunctionSpace_c.f"
#include "atlas_module_PathName_c.f"
#include "atlas_module_Error_c.f"
#include "atlas_module_GridDistribution_c.f"
#include "atlas_module_State_c.f"
#include "atlas_module_Trans_c.f"
#include "atlas_module_Value_c.f"

! -----------------------------------------------------------------------------

function eckit_version()
  character(len=5) :: eckit_version
  eckit_version = c_to_f_string_cptr(atlas__eckit_version())
end function eckit_version

function eckit_git_sha1()
  character(len=40) :: eckit_git_sha1
  eckit_git_sha1 = c_to_f_string_cptr(atlas__eckit_git_sha1())
end function eckit_git_sha1

function eckit_git_sha1_abbrev(length)
  character(len=40) :: eckit_git_sha1_abbrev
  integer(c_int), optional :: length
  integer(c_int) :: opt_length
  opt_length = 7
  if( present(length) ) opt_length = length
  eckit_git_sha1_abbrev = c_to_f_string_cptr(atlas__eckit_git_sha1_abbrev(opt_length))
end function eckit_git_sha1_abbrev

function atlas_version()
  character(len=5) :: atlas_version
  atlas_version = c_to_f_string_cptr(atlas__atlas_version())
end function atlas_version

function atlas_git_sha1()
  character(len=40) :: atlas_git_sha1
  atlas_git_sha1 = c_to_f_string_cptr(atlas__atlas_git_sha1())
end function atlas_git_sha1

function atlas_git_sha1_abbrev(length)
  character(len=40) :: atlas_git_sha1_abbrev
  integer(c_int), optional :: length
  integer(c_int) :: opt_length
  opt_length = 7
  if( present(length) ) opt_length = length
  atlas_git_sha1_abbrev = c_to_f_string_cptr(atlas__atlas_git_sha1_abbrev(opt_length))
end function atlas_git_sha1_abbrev

function atlas_run_name()
  character(len=128) :: atlas_run_name
  atlas_run_name = c_to_f_string_cptr(atlas__run_name())
end function atlas_run_name

function atlas_display_name()
  character(len=128) :: atlas_display_name
  atlas_display_name = c_to_f_string_cptr(atlas__display_name())
end function atlas_display_name

function atlas_rundir()
  character(len=128) :: atlas_rundir
  atlas_rundir = c_to_f_string_cptr(atlas__rundir())
end function atlas_rundir

function atlas_workdir()
  character(len=128) :: atlas_workdir
  atlas_workdir = c_to_f_string_cptr(atlas__workdir())
end function atlas_workdir

function atlas_read_gmsh(filename) result(mesh)
  character(len=*), intent(in) :: filename
  type(atlas_Mesh) :: mesh
  mesh%cpp_object_ptr = atlas__read_gmsh(c_str(filename))
end function atlas_read_gmsh

subroutine atlas_write_gmsh(mesh,filename)
  type(atlas_Mesh), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_gmsh_mesh(mesh%cpp_object_ptr,c_str(filename))
end subroutine atlas_write_gmsh

subroutine atlas_write_gmsh_field(field,function_space,filename,mode)
  type(atlas_Field), intent(in) :: field
  type(atlas_NodesFunctionSpace), intent(in) :: function_space
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_field(field%cpp_object_ptr,function_space%cpp_object_ptr,c_str(filename),mode)
  else
    call atlas__write_gmsh_field(field%cpp_object_ptr,function_space%cpp_object_ptr,c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_field

subroutine atlas_write_gmsh_fieldset(fieldset,function_space,filename,mode)
  type(atlas_FieldSet), intent(in) :: fieldset
  type(atlas_NodesFunctionSpace), intent(in) :: function_space
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_fieldset(fieldset%cpp_object_ptr,function_space%cpp_object_ptr,c_str(filename),mode)
  else
    call atlas__write_gmsh_fieldset(fieldset%cpp_object_ptr,function_space%cpp_object_ptr,c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_fieldset

subroutine atlas_build_parallel_fields(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_parallel_fields(mesh%cpp_object_ptr)
end subroutine atlas_build_parallel_fields

subroutine atlas_build_nodes_parallel_fields(nodes)
  type(atlas_Nodes), intent(inout) :: nodes
  call atlas__build_nodes_parallel_fields(nodes%cpp_object_ptr)
end subroutine atlas_build_nodes_parallel_fields

subroutine atlas_renumber_nodes_glb_idx(nodes)
  type(atlas_Nodes), intent(inout) :: nodes
  call atlas__renumber_nodes_glb_idx(nodes%cpp_object_ptr)
end subroutine atlas_renumber_nodes_glb_idx

subroutine atlas_build_edges_parallel_fields(edges, nodes)
  type(atlas_FunctionSpace), intent(inout) :: edges
  type(atlas_Nodes), intent(inout) :: nodes
  call atlas__build_edges_parallel_fields(edges%cpp_object_ptr,nodes%cpp_object_ptr)
end subroutine atlas_build_edges_parallel_fields

subroutine atlas_build_periodic_boundaries(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_periodic_boundaries(mesh%cpp_object_ptr)
end subroutine atlas_build_periodic_boundaries

subroutine atlas_build_halo(mesh,nelems)
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in) :: nelems
  call atlas__build_halo(mesh%cpp_object_ptr,nelems)
end subroutine atlas_build_halo

subroutine atlas_build_edges(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_edges(mesh%cpp_object_ptr)
end subroutine atlas_build_edges

subroutine atlas_build_pole_edges(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_pole_edges(mesh%cpp_object_ptr)
end subroutine atlas_build_pole_edges

subroutine atlas_build_node_to_edge_connectivity(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_node_to_edge_connectivity(mesh%cpp_object_ptr)
end subroutine atlas_build_node_to_edge_connectivity

subroutine atlas_build_median_dual_mesh(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_median_dual_mesh(mesh%cpp_object_ptr)
end subroutine atlas_build_median_dual_mesh

subroutine atlas_build_centroid_dual_mesh(mesh)
  type(atlas_Mesh), intent(inout) :: mesh
  call atlas__build_centroid_dual_mesh(mesh%cpp_object_ptr)
end subroutine atlas_build_centroid_dual_mesh

subroutine atlas_write_load_balance_report(mesh,filename)
  type(atlas_Mesh), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_load_balance_report(mesh%cpp_object_ptr,c_str(filename))
end subroutine atlas_write_load_balance_report

function atlas_generate_mesh(grid) result(mesh)
  type(atlas_Mesh) :: mesh
  type(atlas_ReducedGrid) :: grid
  mesh%cpp_object_ptr = atlas__generate_mesh(grid%cpp_object_ptr)
end function atlas_generate_mesh

function atlas_generate_mesh_with_distribution(grid,distribution) result(mesh)
  type(atlas_Mesh) :: mesh
  type(atlas_ReducedGrid) :: grid
  type(atlas_GridDistribution) :: distribution
  mesh%cpp_object_ptr = atlas__generate_mesh_with_distribution(grid%cpp_object_ptr,distribution%cpp_object_ptr)
end function atlas_generate_mesh_with_distribution

! -----------------------------------------------------------------------------

end module atlas_module
