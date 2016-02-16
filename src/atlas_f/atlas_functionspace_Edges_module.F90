
module atlas_functionspace_Edges_module

use iso_c_binding, only : c_ptr
use atlas_c_interop, only : c_str, c_to_f_string_cptr, atlas_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Field_module, only: atlas_Field
use atlas_Mesh_module, only: atlas_Mesh
use atlas_mesh_HybridElements_module, only: atlas_mesh_Edges
use atlas_GatherScatter_module, only: atlas_GatherScatter
use atlas_HaloExchange_module, only: atlas_HaloExchange
use atlas_Checksum_module, only: atlas_Checksum

implicit none

private :: c_ptr
private :: c_str, c_to_f_string_cptr, atlas_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_mesh_Edges
private :: atlas_GatherScatter
private :: atlas_HaloExchange
private :: atlas_Checksum
private :: atlas_Mesh

public :: atlas_functionspace_Edges

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_Edges

! Purpose :
! -------
!   *atlas_functionspace_Edges* : Interpretes fields defined in edges

! Methods :
! -------

! Author :
! ------
!   February-2016 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains


  procedure, public :: nb_edges
  procedure, public :: mesh
  procedure, public :: edges

  procedure, private :: create_field_name_kind
  procedure, private :: create_field_name_kind_lev
  procedure, private :: create_field_name_kind_vars
  procedure, private :: create_field_name_kind_lev_vars
  procedure, private :: create_field_name_template
  generic, public :: create_field => &
    & create_field_name_kind, &
    & create_field_name_kind_lev, &
    & create_field_name_kind_vars, &
    & create_field_name_kind_lev_vars, &
    & create_field_name_template

  procedure, private :: create_glb_field_name_kind
  procedure, private :: create_glb_field_name_kind_lev
  procedure, private :: create_glb_field_name_kind_vars
  procedure, private :: create_glb_field_name_kind_lev_vars
  procedure, private :: create_glb_field_name_template
  generic, public :: create_global_field => &
    & create_glb_field_name_kind, &
    & create_glb_field_name_kind_lev, &
    & create_glb_field_name_kind_vars, &
    & create_glb_field_name_kind_lev_vars, &
    & create_glb_field_name_template

  !procedure, private :: halo_exchange_fieldset
  procedure, private :: halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_field !,halo_exchange_fieldset
  procedure, public :: get_halo_exchange

  !procedure, private :: gather_fieldset
  procedure, private :: gather_field
  generic, public :: gather => gather_field !,gather_fieldset
  procedure, public :: get_gather

  !procedure, private :: scatter_fieldset
  procedure, private :: scatter_field
  generic, public :: scatter => scatter_field !,  scatter_fieldset
  procedure, public :: get_scatter

  !procedure, private :: checksum_fieldset
  procedure, private :: checksum_field
  generic, public :: checksum => checksum_field !,checksum_fieldset
  procedure, public :: get_checksum

#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_Edges__final
#endif


END TYPE atlas_functionspace_Edges

interface atlas_functionspace_Edges
  module procedure atlas_functionspace_Edges__cptr
  module procedure atlas_functionspace_Edges__mesh_halo
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

function atlas_functionspace_Edges__cptr(cptr) result(functionspace)
  type(atlas_functionspace_Edges) :: functionspace
  type(c_ptr), intent(in) :: cptr
  call functionspace%reset_c_ptr( cptr )
end function

function atlas_functionspace_Edges__mesh_halo(mesh,halo) result(function_space)
  use atlas_functionspace_Edges_c_binding
  type(atlas_functionspace_Edges) :: function_space
  type(atlas_Mesh), intent(inout) :: mesh
  integer, intent(in), optional :: halo
  if( present(halo) ) then
    function_space = atlas_functionspace_Edges__cptr( &
      & atlas__functionspace__Edges__new(mesh%c_ptr(),halo) )
  else
    function_space = atlas_functionspace_Edges__cptr( &
      & atlas__functionspace__Edges__new_mesh(mesh%c_ptr()) )
  endif
  call function_space%return()
end function


#ifdef FORTRAN_SUPPORTS_FINAL
subroutine atlas_functionspace_Edges__final(this)
  type(atlas_functionspace_Edges), intent(inout) :: this
  call this%final()
end subroutine
#endif

function nb_edges(this)
  use atlas_functionspace_Edges_c_binding
  integer :: nb_edges
  class(atlas_functionspace_Edges), intent(in) :: this
  nb_edges = atlas__functionspace__Edges__nb_edges(this%c_ptr())
end function

function mesh(this)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Mesh) :: mesh
  class(atlas_functionspace_Edges), intent(in) :: this
  call mesh%reset_c_ptr( atlas__functionspace__Edges__mesh(this%c_ptr()) )
end function

function edges(this)
  use atlas_functionspace_Edges_c_binding
  type(atlas_mesh_Edges) :: edges
  class(atlas_functionspace_Edges), intent(in) :: this
  call edges%reset_c_ptr( atlas__functionspace__Edges__edges(this%c_ptr()) )
end function

function create_field_name_kind(this,name,kind) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  field = atlas_Field( atlas__functionspace__Edges__create_field(this%c_ptr(),c_str(name),kind) )
  call field%return()
end function

function create_field_name_kind_lev(this,name,kind,levels) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  integer, intent(in) :: levels
  field = atlas_Field( atlas__functionspace__Edges__create_field_lev(this%c_ptr(),c_str(name),levels,kind) )
  call field%return()
end function

function create_field_name_kind_vars(this,name,kind,vars) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field = atlas_Field( atlas__functionspace__Edges__create_field_vars(this%c_ptr(),c_str(name),vars,size(vars),fortran_ordering,kind) )
  call field%return()
end function

function create_field_name_kind_lev_vars(this,name,kind,levels,vars) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  integer, intent(in) :: levels
  integer, intent(in) :: vars(:)
  integer, parameter :: fortran_ordering = 1
  field = atlas_Field( atlas__functionspace__Edges__create_field_lev_vars(this%c_ptr(),c_str(name),levels,vars,size(vars),fortran_ordering,kind) )
  call field%return()
end function

function create_field_name_template(this,name,template) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: template
  field = atlas_Field( atlas__functionspace__Edges__create_field_template(this%c_ptr(),c_str(name),template%c_ptr()) )
  call field%return()
end function


function create_glb_field_name_kind_lev(this,name,kind) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: kind
  field = atlas_Field( atlas__functionspace__Edges__create_global_field(this%c_ptr(),c_str(name),kind) )
  call field%return()
end function

function create_glb_field_name_kind(this,name,kind,levels) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: levels
  integer, intent(in) :: kind
  field = atlas_Field( atlas__functionspace__Edges__create_global_field_lev(this%c_ptr(),c_str(name),levels,kind) )
  call field%return()
end function


function create_glb_field_name_kind_vars(this,name,kind,vars) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field = atlas_Field( atlas__functionspace__Edges__create_global_field_vars(this%c_ptr(),c_str(name),vars,size(vars),fortran_ordering,kind) )
  call field%return()
end function

function create_glb_field_name_kind_lev_vars(this,name,kind,levels,vars) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: vars(:)
  integer, intent(in) :: levels
  integer, intent(in) :: kind
  integer, parameter :: fortran_ordering = 1
  field = atlas_Field( atlas__functionspace__Edges__create_global_field_lev_vars(this%c_ptr(),c_str(name),levels,vars,size(vars),fortran_ordering,kind) )
  call field%return()
end function

function create_glb_field_name_template(this,name,template) result(field)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_Edges), intent(in) :: this
  character(len=*), intent(in) :: name
  type(atlas_Field) :: template
  field = atlas_Field( atlas__functionspace__Edges__create_global_field_template(this%c_ptr(),c_str(name),template%c_ptr()) )
  call field%return()
end function

#if 0
subroutine halo_exchange_fieldset(this,fieldset)
  use atlas_functionspace_Edges_c_binding
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_FieldSet), intent(inout) :: fieldset
  call atlas__functionspace__Edges__halo_exchange_fieldset(this%c_ptr(),fieldset%c_ptr())
end subroutine
#endif

subroutine halo_exchange_field(this,field)
  use atlas_functionspace_Edges_c_binding
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_Field), intent(inout) :: field
  call atlas__functionspace__Edges__halo_exchange_field(this%c_ptr(),field%c_ptr())
end subroutine

function get_gather(this) result(gather)
  use atlas_functionspace_Edges_c_binding
  type(atlas_GatherScatter) :: gather
  class(atlas_functionspace_Edges), intent(in) :: this
  call gather%reset_c_ptr( atlas__functionspace__Edges__get_gather(this%c_ptr()) )
end function

function get_scatter(this) result(gather)
  use atlas_functionspace_Edges_c_binding
  type(atlas_GatherScatter) :: gather
  class(atlas_functionspace_Edges), intent(in) :: this
  call gather%reset_c_ptr( atlas__functionspace__Edges__get_scatter(this%c_ptr()) )
end function

#if 0
subroutine gather_fieldset(this,local,global)
  use atlas_functionspace_Edges_c_binding
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: local
  type(atlas_FieldSet), intent(inout) :: global
  call atlas__functionspace__Edges__gather_fieldset(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine
#endif

subroutine gather_field(this,local,global)
  use atlas_functionspace_Edges_c_binding
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_Field), intent(in) :: local
  type(atlas_Field), intent(inout) :: global
  call atlas__functionspace__Edges__gather_field(this%c_ptr(),local%c_ptr(),global%c_ptr())
end subroutine

#if 0
subroutine scatter_fieldset(this,global,local)
  use atlas_functionspace_Edges_c_binding
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: global
  type(atlas_FieldSet), intent(inout) :: local
  call atlas__functionspace__Edges__scatter_fieldset(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine
#endif

subroutine scatter_field(this,global,local)
  use atlas_functionspace_Edges_c_binding
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_Field), intent(in) :: global
  type(atlas_Field), intent(inout) :: local
  call atlas__functionspace__Edges__scatter_field(this%c_ptr(),global%c_ptr(),local%c_ptr())
end subroutine

function get_halo_exchange(this) result(halo_exchange)
  use atlas_functionspace_Edges_c_binding
  type(atlas_HaloExchange) :: halo_exchange
  class(atlas_functionspace_Edges), intent(in) :: this
  call halo_exchange%reset_c_ptr( atlas__functionspace__Edges__get_halo_exchange(this%c_ptr()) )
end function

function get_checksum(this) result(checksum)
  use atlas_functionspace_Edges_c_binding
  type(atlas_Checksum) :: checksum
  class(atlas_functionspace_Edges), intent(in) :: this
  call checksum%reset_c_ptr( atlas__functionspace__Edges__get_checksum(this%c_ptr()) )
end function

#if 0
function checksum_fieldset(this,fieldset) result(checksum)
  use atlas_functionspace_Edges_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__functionspace__Edges__checksum_fieldset(this%c_ptr(),fieldset%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function
#endif

function checksum_field(this,field) result(checksum)
  use atlas_functionspace_Edges_c_binding
  character(len=:), allocatable :: checksum
  class(atlas_functionspace_Edges), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(c_ptr) :: checksum_cptr
  integer :: checksum_size, checksum_allocated
  call atlas__functionspace__Edges__checksum_field(this%c_ptr(),field%c_ptr(),checksum_cptr,checksum_size,checksum_allocated)
  allocate(character(len=checksum_size) :: checksum )
  checksum = c_to_f_string_cptr(checksum_cptr)
  if( checksum_allocated == 1 ) call atlas_free(checksum_cptr)
end function


end module atlas_functionspace_Edges_module

