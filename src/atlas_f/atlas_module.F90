! (C) Copyright 2013-2015 ECMWF.

#include "atlas/atlas_f.h"

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

use atlas_mpi_module

use atlas_FunctionSpace_module, only: &
    & atlas_FunctionSpace
use atlas_Field_module, only: &
    & atlas_Field, &
    & atlas_real, &
    & atlas_integer, &
    & atlas_logical
use atlas_FieldSet_module, only: &
    & atlas_FieldSet
use atlas_State_module, only: &
    & atlas_State
use atlas_Config_module, only: &
    & atlas_config
use atlas_JSON_module, only: &
    & atlas_JSON, &
    & atlas_PathName
use atlas_Metadata_module, only: &
    & atlas_Metadata
use atlas_Logging_module, only: &
    & atlas_log, &
    & atlas_Logger, &
    & atlas_LogChannel, &
    & ATLAS_LOG_CATEGORY_ERROR, &
    & ATLAS_LOG_CATEGORY_WARNING, &
    & ATLAS_LOG_CATEGORY_INFO, &
    & ATLAS_LOG_CATEGORY_DEBUG, &
    & ATLAS_LOG_CATEGORY_STATS
use atlas_Error_module, only: &
    & atlas_CodeLocation, &
    & atlas_code_location_str, &
    & atlas_code_location, &
    & atlas_abort, &
    & atlas_throw_exception, &
    & atlas_throw_notimplemented, &
    & atlas_throw_outofrange, &
    & atlas_throw_seriousbug, &
    & atlas_throw_usererror, &
    & atlas_throw_assertionfailed, &
    & atlas_err, &
    & atlas_noerr, &
    & atlas_err_clear, &
    & atlas_err_success, &
    & atlas_err_code, &
    & atlas_err_msg, &
    & atlas_err_set_aborts, &
    & atlas_err_set_throws, &
    & atlas_err_set_backtrace, &
    & atlas_err_cleared, &
    & atlas_err_noerr, &
    & atlas_err_exception, &
    & atlas_err_usererror, &
    & atlas_err_seriousbug, &
    & atlas_err_notimplemented, &
    & atlas_err_assertionfailed, &
    & atlas_err_badparameter, &
    & atlas_err_outofrange, &
    & atlas_err_stop, &
    & atlas_err_abort, &
    & atlas_err_cancel, &
    & atlas_err_readerror, &
    & atlas_err_writeerror, &
    & atlas_err_unknown
use atlas_HybridElements_module, only: &
    & atlas_HybridElements
use atlas_mesh_Edges_module, only: &
    & atlas_mesh_Edges
use atlas_mesh_Cells_module, only: &
    & atlas_mesh_Cells
use atlas_Elements_module, only: &
    & atlas_Elements
use atlas_ElementType_module, only: &
    & atlas_ElementType, &
    & atlas_Triangle, &
    & atlas_Quadrilateral, &
    & atlas_Line
use atlas_Connectivity_module, only: &
    & atlas_Connectivity, &
    & atlas_MultiBlockConnectivity, &
    & atlas_BlockConnectivity
use atlas_mesh_Nodes_module, only: &
    & atlas_mesh_Nodes
use atlas_HaloExchange_module, only: &
    & atlas_HaloExchange
use atlas_GatherScatter_module, only: &
    & atlas_GatherScatter
use atlas_Checksum_module, only: &
    & atlas_Checksum
use atlas_Mesh_module, only: &
    & atlas_Mesh
use atlas_Grid_module, only: &
    & atlas_Grid, &
    & atlas_grid_Structured, &
    & atlas_grid_CustomStructured, &
    & atlas_grid_ReducedGaussian, &
    & atlas_grid_RegularGaussian, &
    & atlas_grid_RegularLonLat, &
    & atlas_grid_ShiftedLonLat, &
    & atlas_grid_ShiftedLon, &
    & atlas_grid_ShiftedLat
use atlas_functionspace_EdgeColumns_module, only: &
    & atlas_functionspace_EdgeColumns
use atlas_functionspace_NodeColumns_module, only: &
    & atlas_functionspace_NodeColumns
use atlas_functionspace_StructuredColumns_module, only: &
    & atlas_functionspace_StructuredColumns
use atlas_functionspace_Spectral_module, only: &
    & atlas_functionspace_Spectral
use atlas_Trans_module, only : &
    & atlas_Trans
use atlas_kinds_module, only: &
    & ATLAS_KIND_GIDX, &
    & ATLAS_KIND_IDX
use atlas_GridDistribution_module, only: &
    & atlas_GridDistribution
use atlas_MeshGenerator_module, only: &
    & atlas_MeshGenerator, &
    & atlas_meshgenerator_Structured
use atlas_Method_module, only: &
    & atlas_Method
use atlas_fvm_module, only: &
    & atlas_fvm_Method
use atlas_Nabla_module, only: &
    & atlas_Nabla
use atlas_resource_module, only: &
    & atlas_resource, &
    & atlas_resource_set
use atlas_Value_module, only: &
    & atlas_Value
use atlas_mesh_actions_module, only: &
    & atlas_build_parallel_fields, &
    & atlas_build_nodes_parallel_fields, &
    & atlas_build_edges_parallel_fields, &
    & atlas_build_periodic_boundaries, &
    & atlas_build_halo, &
    & atlas_build_edges, &
    & atlas_build_pole_edges, &
    & atlas_build_node_to_edge_connectivity, &
    & atlas_build_median_dual_mesh, &
    & atlas_write_load_balance_report, &
    & atlas_renumber_nodes_glb_idx
use atlas_io_module, only: &
    & atlas_read_gmsh, &
    & atlas_write_gmsh, &
    & atlas_write_gmsh_field, &
    & atlas_write_gmsh_fieldset, &
    & atlas_Output, &
    & atlas_output_Gmsh

implicit none

public

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

! =============================================================================
CONTAINS
! =============================================================================


subroutine atlas_init( mpi_comm )
  use, intrinsic :: iso_c_binding, only : c_ptr
  use atlas_atlas_c_binding
  use atlas_mpi_module, only :  atlas_mpi_comm_attach_fortran_communicator
  use atlas_c_interop
  integer, save :: argc
  type(c_ptr), save :: argv(15)
  integer, intent(in), optional :: mpi_comm
  call get_c_arguments(argc,argv)
  if( present(mpi_comm) ) then
    call atlas_mpi_comm_attach_fortran_communicator(mpi_comm)
  endif
  call atlas__atlas_init(argc,argv)
  atlas_log = atlas_Logger()
end subroutine

subroutine atlas_finalize()
  use atlas_c_interop
  use atlas_atlas_c_binding
  call atlas__atlas_finalize()
end subroutine

! -----------------------------------------------------------------------------

function eckit_version()
  use atlas_atlas_c_binding
  use atlas_c_interop
  character(len=5) :: eckit_version
  eckit_version = c_to_f_string_cptr(atlas__eckit_version())
end function eckit_version

function eckit_git_sha1()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=40) :: eckit_git_sha1
  eckit_git_sha1 = c_to_f_string_cptr(atlas__eckit_git_sha1())
end function eckit_git_sha1

function eckit_git_sha1_abbrev(length)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_atlas_c_binding
  use atlas_c_interop
  character(len=40) :: eckit_git_sha1_abbrev
  integer(c_int), optional :: length
  integer(c_int) :: opt_length
  opt_length = 7
  if( present(length) ) opt_length = length
  eckit_git_sha1_abbrev = c_to_f_string_cptr(atlas__eckit_git_sha1_abbrev(opt_length))
end function eckit_git_sha1_abbrev

function atlas_version()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=5) :: atlas_version
  atlas_version = c_to_f_string_cptr(atlas__atlas_version())
end function atlas_version

function atlas_git_sha1()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=40) :: atlas_git_sha1
  atlas_git_sha1 = c_to_f_string_cptr(atlas__atlas_git_sha1())
end function atlas_git_sha1

function atlas_git_sha1_abbrev(length)
  use atlas_c_interop
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_atlas_c_binding
  character(len=40) :: atlas_git_sha1_abbrev
  integer(c_int), optional :: length
  integer(c_int) :: opt_length
  opt_length = 7
  if( present(length) ) opt_length = length
  atlas_git_sha1_abbrev = c_to_f_string_cptr(atlas__atlas_git_sha1_abbrev(opt_length))
end function atlas_git_sha1_abbrev

function atlas_run_name()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=128) :: atlas_run_name
  atlas_run_name = c_to_f_string_cptr(atlas__run_name())
end function atlas_run_name

function atlas_display_name()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=128) :: atlas_display_name
  atlas_display_name = c_to_f_string_cptr(atlas__display_name())
end function atlas_display_name

function atlas_rundir()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=128) :: atlas_rundir
  atlas_rundir = c_to_f_string_cptr(atlas__rundir())
end function atlas_rundir

function atlas_workdir()
  use atlas_c_interop
  use atlas_atlas_c_binding
  character(len=128) :: atlas_workdir
  atlas_workdir = c_to_f_string_cptr(atlas__workdir())
end function atlas_workdir


! -----------------------------------------------------------------------------

end module atlas_module
