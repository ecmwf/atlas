! (C) Copyright 2013-2015 ECMWF.

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
use atlas_object_module
use atlas_refcounted_module, only: &
    & atlas_refcounted, &
    & atlas_refcounted_fortran
use atlas_FunctionSpace_module, only: &
    & atlas_FunctionSpace
use atlas_Field_module, only: &
    & atlas_Field, &
    & atlas_real, &
    & atlas_integer, &
    & atlas_logical
use atlas_FieldSet_module, only: &
    & atlas_FieldSet
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
use atlas_mesh_HybridElements_module, only: &
    & atlas_mesh_HybridElements, &
    & atlas_mesh_Cells, &
    & atlas_mesh_Edges
use atlas_mesh_Elements_module, only: &
    & atlas_mesh_Elements
use atlas_mesh_ElementType_module, only: &
    & atlas_mesh_ElementType, &
    & atlas_mesh_Triangle, &
    & atlas_mesh_Quadrilateral, &
    & atlas_mesh_Line
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
    & atlas_ReducedGrid, &
    & atlas_ReducedGaussianGrid, &
    & atlas_GaussianGrid, &
    & atlas_LonLatGrid
use atlas_functionspace_Edges_module, only: &
    & atlas_functionspace_Edges
use atlas_functionspace_Nodes_module, only: &
    & atlas_functionspace_Nodes
use atlas_kinds_module, only: &
    & ATLAS_KIND_GIDX, &
    & ATLAS_KIND_IDX
use atlas_GridDistribution_module, only: &
    & atlas_GridDistribution
use atlas_MeshGenerator_module, only: &
    & atlas_MeshGenerator, &
    & atlas_ReducedGridMeshGenerator
use atlas_numerics_Method_module, only: &
    & atlas_numerics_Method
use atlas_numerics_fvm_module, only: &
    & atlas_numerics_fvm_Method
use atlas_numerics_Nabla_module, only: &
    & atlas_numerics_Nabla


use atlas_actions_module

#if !DEPRECATE_OLD_FUNCTIONSPACE
use atlas_deprecated_functionspace_module, only: atlas_deprecated_FunctionSpace
#endif

use atlas_atlas_c_binding
implicit none

private :: atlas_object
private :: atlas_refcounted
private :: atlas_refcounted_fortran
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

public


#include "atlas_module_functionspace_ReducedGridPoint_i.f"
#include "atlas_module_functionspace_Spectral_i.f"
#include "atlas_module_State_i.f"
#include "atlas_module_Trans_i.f"
#include "atlas_module_Value_i.f"

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

#include "atlas_module_functionspace_ReducedGridPoint_c.f"
#include "atlas_module_functionspace_Spectral_c.f"
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


! -----------------------------------------------------------------------------

end module atlas_module
