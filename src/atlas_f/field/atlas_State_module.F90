
module atlas_State_module

use fckit_owned_object_module, only: fckit_owned_object
use atlas_Field_module, only: atlas_Field

implicit none

private :: fckit_owned_object
private :: atlas_Field

public :: atlas_State

private

!----------------------------!
! atlas_State                !
!----------------------------!

! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_State

! Purpose :
! -------
!   *atlas_State* :
!       Container type of fields that are defined on the same points
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done
!       Describes interpolation between nodes

! Methods :
! -------
!   name : The name or tag this function space was created with
!   create_field : Create a new real field in this function space with given name
!   remove : Remove a field with given name
!   field : Access to a field with given name
!   parallelise : Setup halo-exchange information
!   halo_exchange : Perform halo exchange on field_data

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

!-- Field
  procedure, public :: add    => atlas_State__add
  procedure, public :: remove => atlas_State__remove
  procedure, public :: has    => atlas_State__has
  procedure, public :: size    => atlas_State__size
  procedure, private :: field_by_name  => atlas_State__field_by_name
  procedure, private :: field_by_index => atlas_State__field_by_index
  generic, public :: field => field_by_name, field_by_index
  procedure, public :: metadata => atlas_State__metadata
END TYPE atlas_State

interface atlas_State
  module procedure atlas_State__new
  module procedure atlas_State__generate
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! State routines


function atlas_State__new() result(this)
  use atlas_state_c_binding
  type(atlas_State) :: this
  call this%reset_c_ptr( atlas__State__new() )
  call this%return()
end function

function atlas_State__generate(generator, params) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_state_c_binding
  use atlas_Config_module, only: atlas_Config
  type(atlas_State) :: this
  character(len=*), intent(in) :: generator
  class(atlas_Config), intent(in), optional :: params

  type(atlas_Config) :: p

  call this%reset_c_ptr( atlas__State__new() )

  if( present(params) ) then
    call atlas__State__initialize(this%c_ptr(),c_str(generator),params%c_ptr())
  else
    p = atlas_Config()
    call atlas__State__initialize(this%c_ptr(),c_str(generator),p%c_ptr())
    call p%final()
  endif
  call this%return()
end function

subroutine atlas_State__add(this,field)
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  class(atlas_Field), intent(in) :: field
  call atlas__State__add(this%c_ptr(),field%c_ptr())
end subroutine
 
subroutine atlas_State__remove(this,name)
  use fckit_c_interop_module, only: c_str
  use atlas_state_c_binding
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  call atlas__State__remove(this%c_ptr(),c_str(name))
end subroutine

function atlas_State__has(this,name) result(has)
  use fckit_c_interop_module, only: c_str
  use atlas_state_c_binding
  logical :: has
  class(atlas_State), intent(in) :: this
  character(len=*), intent(in) :: name
  integer :: has_int
  has_int = atlas__State__has(this%c_ptr(),c_str(name))
  has = .False.
  if( has_int == 1 ) has = .True.
end function

function atlas_State__size(this) result(size)
  use atlas_state_c_binding
  integer :: size
  class(atlas_State), intent(in) :: this
  size = atlas__State__size(this%c_ptr())
end function

function atlas_State__field_by_name(this,name) result(field)
  use fckit_c_interop_module, only: c_str
  use atlas_state_c_binding
  type(atlas_Field) :: field
  class(atlas_State), intent(inout) :: this
  character(len=*), intent(in) :: name
  field = atlas_Field( atlas__State__field_by_name(this%c_ptr(),c_str(name)) )
  call field%return()
end function

function atlas_State__field_by_index(this,index) result(field)
  use atlas_state_c_binding
  type(atlas_Field) :: field
  class(atlas_State), intent(in) :: this
  integer, intent(in) :: index
  field = atlas_Field( atlas__State__field_by_index(this%c_ptr(),index-1) )
  call field%return()
end function

function atlas_State__metadata(this) result(metadata)
  use atlas_state_c_binding
  use atlas_Metadata_module, only: atlas_Metadata
  type(atlas_Metadata) :: metadata
  class(atlas_State), intent(in) :: this
  call metadata%reset_c_ptr( atlas__State__metadata(this%c_ptr()) )
end function

! ----------------------------------------------------------------------------------------

end module atlas_State_module
