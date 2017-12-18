! (C) Copyright 2013-2015 ECMWF.

#include "atlas/atlas_f.h"

module atlas_module

! Purpose :
! -------
!    *atlas* : Object oriented library for flexible parallel data structures
!              for structured grids and unstructured meshes
!
!------------------------------------------------------------------------------

use atlas_mpi_module

use atlas_Field_module, only: &
    & atlas_Field, &
    & atlas_real, &
    & atlas_integer, &
    & atlas_logical
use atlas_FunctionSpace_module, only: &
    & atlas_FunctionSpace
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
use atlas_Interpolation_module, only: &
    & atlas_Interpolation
use atlas_GatherScatter_module, only: &
    & atlas_GatherScatter
use atlas_Checksum_module, only: &
    & atlas_Checksum
use atlas_Mesh_module, only: &
    & atlas_Mesh
use atlas_Grid_module, only: &
    & atlas_Grid , &
    & atlas_StructuredGrid, &
    & atlas_GaussianGrid, &
    & atlas_ReducedGaussianGrid, &
    & atlas_RegularGaussianGrid, &
    & atlas_RegularLonLatGrid
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
    & ATLAS_KIND_IDX, &
    & ATLAS_KIND_REAL64, &
    & ATLAS_KIND_REAL32, &
    & ATLAS_KIND_INT64, &
    & ATLAS_KIND_INT32
use atlas_GridDistribution_module, only: &
    & atlas_GridDistribution
use atlas_Partitioner_module, only: &
    & atlas_Partitioner, &
    & atlas_MatchingMeshPartitioner
use atlas_MeshGenerator_module, only: &
    & atlas_MeshGenerator
use atlas_Method_module, only: &
    & atlas_Method
use atlas_fvm_module, only: &
    & atlas_fvm_Method
use atlas_Nabla_module, only: &
    & atlas_Nabla
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
use atlas_output_module, only: &
    & atlas_Output, &
    & atlas_output_Gmsh

use fckit_log_module,  only: atlas_log => fckit_log

implicit none

public

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

type, private :: atlas_Library_type
contains
  procedure, public, nopass :: initialise => atlas_init
  procedure, public, nopass :: finalise   => atlas_finalise
  procedure, public, nopass :: version    => atlas_version
  procedure, public, nopass :: gitsha1    => atlas_git_sha1_abbrev
end type

type(atlas_library_type), public :: atlas_library

  type, private :: eckit_Library_type
  contains
    procedure, public, nopass :: version    => eckit_version
    procedure, public, nopass :: gitsha1    => eckit_git_sha1_abbrev
  end type

type(eckit_library_type), public :: eckit_library

! =============================================================================
CONTAINS
! =============================================================================

! -----------------------------------------------------------------------------

subroutine atlas_init( comm )
  use atlas_library_c_binding
  use iso_fortran_env, only : stdout => output_unit
  use fckit_main_module, only: fckit_main
  use fckit_mpi_module, only : fckit_mpi_comm
  use atlas_mpi_module, only : atlas_mpi_set_comm

  integer, intent(in), optional :: comm
  type(fckit_mpi_comm) :: world
  integer :: output_unit

  if( .not. fckit_main%ready() ) then
    call fckit_main%initialise()

    if( fckit_main%taskID() == 0 ) then
      call atlas_log%set_fortran_unit(stdout,style=atlas_log%PREFIX)
    else
      call atlas_log%reset()
    endif

    call atlas_log%debug("--> Only MPI rank 0 is logging. Please initialise fckit_main"//&
      &                  "    before to avoid this behaviour",flush=.true.);

  endif

  if( present(comm) ) then
    call atlas_mpi_set_comm(comm)
  endif
  call atlas__atlas_init_noargs()

end subroutine



subroutine atlas_finalise()
  use fckit_c_interop_module
  use atlas_library_c_binding
  call atlas__atlas_finalize()
end subroutine

function eckit_version()
  use atlas_library_c_binding
  use fckit_c_interop_module
  character(len=40) :: eckit_version
  eckit_version = c_ptr_to_string(atlas__eckit_version())
end function eckit_version

function eckit_git_sha1()
  use fckit_c_interop_module
  use atlas_library_c_binding
  character(len=40) :: eckit_git_sha1
  eckit_git_sha1 = c_ptr_to_string(atlas__eckit_git_sha1())
end function eckit_git_sha1

function eckit_git_sha1_abbrev(length)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_library_c_binding
  use fckit_c_interop_module
  character(len=40) :: eckit_git_sha1_abbrev
  integer(c_int), optional :: length
  integer(c_int) :: opt_length
  opt_length = 7
  if( present(length) ) opt_length = length
  eckit_git_sha1_abbrev = c_ptr_to_string(atlas__eckit_git_sha1_abbrev(opt_length))
end function eckit_git_sha1_abbrev

function atlas_version()
  use fckit_c_interop_module
  use atlas_library_c_binding
  character(len=40) :: atlas_version
  atlas_version = c_ptr_to_string(atlas__atlas_version())
end function atlas_version

function atlas_git_sha1()
  use fckit_c_interop_module
  use atlas_library_c_binding
  character(len=40) :: atlas_git_sha1
  atlas_git_sha1 = c_ptr_to_string(atlas__atlas_git_sha1())
end function atlas_git_sha1

function atlas_git_sha1_abbrev(length)
  use fckit_c_interop_module
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_library_c_binding
  character(len=40) :: atlas_git_sha1_abbrev
  integer(c_int), optional :: length
  integer(c_int) :: opt_length
  opt_length = 7
  if( present(length) ) opt_length = length
  atlas_git_sha1_abbrev = c_ptr_to_string(atlas__atlas_git_sha1_abbrev(opt_length))
end function atlas_git_sha1_abbrev

! -----------------------------------------------------------------------------

end module atlas_module
