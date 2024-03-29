! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_functionspace_CellColumns_module

use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Mesh_module, only: atlas_Mesh
use atlas_mesh_Cells_module, only: atlas_mesh_Cells
use atlas_GatherScatter_module, only: atlas_GatherScatter
use atlas_HaloExchange_module, only: atlas_HaloExchange
use atlas_Checksum_module, only: atlas_Checksum
use atlas_Config_module, only: atlas_Config
use atlas_kinds_module, only: ATLAS_KIND_GIDX

implicit none

private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_mesh_Cells
private :: atlas_GatherScatter
private :: atlas_HaloExchange
private :: atlas_Checksum
private :: atlas_Mesh
private :: atlas_Config
private :: ATLAS_KIND_GIDX

public :: atlas_functionspace_CellColumns

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_CellColumns

! Purpose :
! -------
!   *atlas_functionspace_CellColumns* : Interpretes fields defined in cells

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains


  procedure, public :: nb_cells
  procedure, public :: mesh
  procedure, public :: cells

  procedure, public :: get_halo_exchange

  procedure, private :: gather_fieldset
  procedure, private :: gather_field
  generic, public :: gather => gather_fieldset, gather_field
  procedure, public :: get_gather

  procedure, private :: scatter_fieldset
  procedure, private :: scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field
  procedure, public :: get_scatter

  procedure, private :: checksum_fieldset
  procedure, private :: checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field
  procedure, public :: get_checksum

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_functionspace_CellColumns__final_auto
#endif

END TYPE atlas_functionspace_CellColumns

interface atlas_functionspace_CellColumns
  module procedure constructor__cptr
  module procedure constructor
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

!------------------------------------------------------------------------------

function constructor__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_functionspace_CellColumns) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

!------------------------------------------------------------------------------

function constructor(mesh,halo,levels) result(this)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_functionspace_CellColumns) :: this
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  integer, intent(in), optional :: levels
  type(atlas_Config) :: config
  config = atlas_Config()
  if( present(halo) )   call config%set("halo",halo)
  if( present(levels) ) call config%set("levels",levels)
  call this%reset_c_ptr( atlas__CellsFunctionSpace__new(mesh%CPTR_PGIBUG_A,config%CPTR_PGIBUG_B) )
  call config%final()
  call this%return()
end function

!------------------------------------------------------------------------------

function nb_cells(this)
  use atlas_functionspace_CellColumns_c_binding
  integer :: nb_cells
  class(atlas_functionspace_CellColumns), intent(in) :: this
  nb_cells = atlas__CellsFunctionSpace__nb_cells(this%CPTR_PGIBUG_A)
end function

!------------------------------------------------------------------------------

function mesh(this)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_functionspace_CellColumns), intent(in) :: this
  call mesh%reset_c_ptr( atlas__CellsFunctionSpace__mesh(this%CPTR_PGIBUG_A) )
  call mesh%return()
end function

!------------------------------------------------------------------------------

function cells(this)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_mesh_Cells) :: cells
  class(atlas_functionspace_CellColumns), intent(in) :: this
  call cells%reset_c_ptr( atlas__CellsFunctionSpace__cells(this%CPTR_PGIBUG_A) )
  call cells%return()
end function

!------------------------------------------------------------------------------

function get_gather(this) result(gather)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_GatherScatter) :: gather
  class(atlas_functionspace_CellColumns), intent(in) :: this
  call gather%reset_c_ptr( atlas__CellsFunctioNSpace__get_gather(this%CPTR_PGIBUG_A) )
end function

!------------------------------------------------------------------------------

function get_scatter(this) result(gather)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_GatherScatter) :: gather
  class(atlas_functionspace_CellColumns), intent(in) :: this
  call gather%reset_c_ptr( atlas__CellsFunctioNSpace__get_scatter(this%CPTR_PGIBUG_A) )
end function

!------------------------------------------------------------------------------

subroutine gather_fieldset(this,local,global)
  use atlas_functionspace_CellColumns_c_binding
  class(atlas_functionspace_CellColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__CellsFunctionSpace__gather_fieldset(this%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

subroutine gather_field(this,local,global)
  use atlas_functionspace_CellColumns_c_binding
  class(atlas_functionspace_CellColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__CellsFunctionSpace__gather_field(this%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

subroutine scatter_fieldset(this,global,local)
  use atlas_functionspace_CellColumns_c_binding
  class(atlas_functionspace_CellColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__CellsFunctionSpace__scatter_fieldset(this%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

subroutine scatter_field(this,global,local)
  use atlas_functionspace_CellColumns_c_binding
  class(atlas_functionspace_CellColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__CellsFunctionSpace__scatter_field(this%CPTR_PGIBUG_A,global%CPTR_PGIBUG_A,local%CPTR_PGIBUG_A)
end subroutine

!------------------------------------------------------------------------------

function get_halo_exchange(this) result(halo_exchange)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_HaloExchange) :: halo_exchange
  class(atlas_functionspace_CellColumns), intent(in) :: this
  call halo_exchange%reset_c_ptr( atlas__CellsFunctioNSpace__get_halo_exchange(this%CPTR_PGIBUG_A) )
end function

!------------------------------------------------------------------------------

function get_checksum(this) result(checksum)
  use atlas_functionspace_CellColumns_c_binding
  type(atlas_Checksum) :: checksum
  class(atlas_functionspace_CellColumns), intent(in) :: this
  call checksum%reset_c_ptr( atlas__CellsFunctioNSpace__get_checksum(this%CPTR_PGIBUG_A) )
end function

!------------------------------------------------------------------------------

function checksum_fieldset(this,fieldset) result(checksum)
  use atlas_functionspace_CellColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_CellColumns), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__CellsFunctionSpace__checksum_fieldset( &
    this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

!------------------------------------------------------------------------------

function checksum_field(this,field) result(checksum)
  use atlas_functionspace_CellColumns_c_binding
  use, intrinsic :: iso_c_binding, only : c_ptr
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_CellColumns), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__CellsFunctionSpace__checksum_field( &
    this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A,checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_ptr_to_string(checksum_cptr)
  if( checksum_allocated == 1 ) call c_ptr_free(checksum_cptr)
end function

!------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_functionspace_CellColumns__final_auto(this)
  type(atlas_functionspace_CellColumns), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_functionspace_CellColumns__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!------------------------------------------------------------------------------

end module atlas_functionspace_CellColumns_module

