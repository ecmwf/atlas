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
use atlas_field_c_binding
use atlas_fieldset_c_binding
use atlas_functionspace_c_binding
use atlas_mesh_c_binding
use atlas_metadata_c_binding
use atlas_haloexchange_c_binding
use atlas_gmsh_c_binding
use atlas_BuildPeriodicBoundaries_c_binding
use atlas_BuildEdges_c_binding
use atlas_BuildDualMesh_c_binding
use atlas_BuildParallelFields_c_binding
use atlas_BuildHalo_c_binding
implicit none

integer, private, parameter :: MAX_STR_LEN = 255
integer, private, parameter :: FIELD_NB_VARS = -1
integer, private, parameter :: KIND_INT32  = -4
integer, private, parameter :: KIND_REAL32 = 4
integer, private, parameter :: KIND_REAL64 = 8
integer, private, parameter :: wp = c_double ! working precision

!------------------------------------------------------------------------------
TYPE :: Mesh_type

! Purpose :
! -------
!   *Mesh* : Container type holding an entire mesh

! Methods :
! -------
!   add_function_space : Add a new function space
!   function_space : Access the function space with given name

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure :: add_function_space => Mesh__add_function_space
  procedure :: function_space => Mesh__function_space
END TYPE Mesh_type

interface new_Mesh
  module procedure new_Mesh
end interface
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
TYPE :: FunctionSpace_type

! Purpose :
! -------
!   *FunctionSpace* : 
!       Container type of fields that are defined on the same points
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done
!       Describes interpolation between nodes

! Methods :
! -------
!   name : The name or tag this function space was created with
!   create_field : Create a new real field in this function space with given name
!   remove_field : Remove a field with given name
!   field : Access to a field with given name
!   parallelise : Setup halo-exchange information
!   halo_exchange : Perform halo exchange on field_data

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure :: name => FunctionSpace__name
  procedure :: dof => FunctionSpace__dof
  procedure :: glb_dof => FunctionSpace__glb_dof
  procedure :: create_field => FunctionSpace__create_field
  procedure :: remove_field => FunctionSpace__remove_field
  procedure :: bounds => FunctionSpace__bounds
  procedure :: field => FunctionSpace__field
  procedure, public :: has_field => FunctionSpace__has_field
  procedure :: parallelise => FunctionSpace__parallelise
  procedure, private :: FunctionSpace__halo_exchange_int32_r1
  procedure, private :: FunctionSpace__halo_exchange_int32_r2
  procedure, private :: FunctionSpace__halo_exchange_int32_r3
  procedure, private :: FunctionSpace__halo_exchange_real32_r1
  procedure, private :: FunctionSpace__halo_exchange_real32_r2
  procedure, private :: FunctionSpace__halo_exchange_real32_r3
  procedure, private :: FunctionSpace__halo_exchange_real64_r1
  procedure, private :: FunctionSpace__halo_exchange_real64_r2
  procedure, private :: FunctionSpace__halo_exchange_real64_r3
  procedure, private :: FunctionSpace__halo_exchange_real64_r4
  procedure :: get_halo_exchange => FunctionSpace__get_halo_exchange
  generic :: halo_exchange => &
      & FunctionSpace__halo_exchange_int32_r1, &
      & FunctionSpace__halo_exchange_int32_r2, &
      & FunctionSpace__halo_exchange_int32_r3, &
      & FunctionSpace__halo_exchange_real32_r1, &
      & FunctionSpace__halo_exchange_real32_r2, &
      & FunctionSpace__halo_exchange_real32_r3, &
      & FunctionSpace__halo_exchange_real64_r1, &
      & FunctionSpace__halo_exchange_real64_r2, &
      & FunctionSpace__halo_exchange_real64_r3, &
      & FunctionSpace__halo_exchange_real64_r4
      procedure, private :: FunctionSpace__gather_real32_r1
      procedure, private :: FunctionSpace__gather_real32_r2
      procedure, private :: FunctionSpace__gather_real32_r3
  procedure, private :: FunctionSpace__gather_real64_r1
  procedure, private :: FunctionSpace__gather_real64_r2
  procedure, private :: FunctionSpace__gather_real64_r3
  procedure, private :: FunctionSpace__gather_int32_r1
  procedure, private :: FunctionSpace__gather_int32_r2
  generic :: gather => &
      & FunctionSpace__gather_real32_r1, &
      & FunctionSpace__gather_real32_r2, &
      & FunctionSpace__gather_real32_r3, &
      & FunctionSpace__gather_real64_r1, &
      & FunctionSpace__gather_real64_r2, &
      & FunctionSpace__gather_real64_r3, &
      & FunctionSpace__gather_int32_r1, &
      & FunctionSpace__gather_int32_r2
END TYPE FunctionSpace_type

interface new_FunctionSpace
  module procedure new_FunctionSpace
end interface
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
TYPE :: Field_type

! Purpose :
! -------
!   *Field* : Object containing field data and Metadata

! Methods :
! -------
!   name : The name or tag this field was created with
!   data : Return the field as a fortran array of specified shape
!   Metadata : Return object that can contain a variety of Metadata

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure :: name => Field__name
  procedure :: data_type => Field__data_type
  procedure :: nb_vars => Field__nb_vars
  procedure :: metadata => Field__metadata
  procedure :: function_space => Field__function_space
  procedure, private :: access_data1_integer => Field__access_data1_integer
  procedure, private :: access_data2_integer => Field__access_data2_integer
  procedure, private :: access_data3_integer => Field__access_data3_integer
  procedure, private :: access_data1_real32 => Field__access_data1_real32
  procedure, private :: access_data2_real32 => Field__access_data2_real32
  procedure, private :: access_data3_real32 => Field__access_data3_real32
  procedure, private :: access_data1_real64 => Field__access_data1_real64
  procedure, private :: access_data2_real64 => Field__access_data2_real64
  procedure, private :: access_data3_real64 => Field__access_data3_real64
  generic :: access_data => access_data1_real32, access_data2_real32, access_data3_real32, &
    & access_data1_real64, access_data2_real64, access_data3_real64, &
    & access_data1_integer, access_data2_integer, access_data3_integer
  procedure :: data1 => Field__data1_wp
  procedure :: data2 => Field__data2_wp
  procedure :: data3 => Field__data3_wp
END TYPE Field_type
!------------------------------------------------------------------------------





!------------------------------------------------------------------------------
TYPE :: FieldSet_type

! Purpose :
! -------
!   *FieldSet* : Object that groups Fields that go together
!       Fields can belong to several fieldsets simultaneously.
!       The actual ownership of the field lies in a FunctionSpace

! Methods :
! -------
!   add_field : The name or tag this field was created with
!   field : Return the field as a fortran array of specified shape
!   get_array : allocate a list of fields contained in the fieldset

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure, public :: size => FieldSet__size
  procedure, public :: add_field => FieldSet__add_field
  procedure, public :: has_field => FieldSet__has_field
  procedure, private :: field_by_name => FieldSet__field_by_name
  procedure, private :: field_by_idx => FieldSet__field_by_idx
  generic :: field => field_by_name, field_by_idx
  procedure, public :: get_array => FieldSet__fields
END TYPE FieldSet_type
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
TYPE :: Metadata_type

! Purpose :
! -------
!   *Metadata* : Container of Metadata, parameters or attributes
!       The Metadata are added as key, value pairs

! Methods :
! -------
!   add : Add a new property with given key and value
!   set : Modify a property with given key and value
!   get : Return a property value for given key

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure, private :: add_logical => Metadata__add_logical
  procedure, private :: add_integer => Metadata__add_integer
  procedure, private :: add_real32 => Metadata__add_real32
  procedure, private :: add_real64 => Metadata__add_real64
  procedure, private :: add_string => Metadata__add_string
  generic :: add => add_logical, add_integer, add_real32, add_real64, add_string
  generic :: set => add_logical, add_integer, add_real32, add_real64, add_string
  procedure :: get_integer => Metadata__get_integer
  procedure :: get_logical => Metadata__get_logical
  procedure :: get_real32 => Metadata__get_real32
  procedure :: get_real64 => Metadata__get_real64
  procedure :: get_string => Metadata__get_string
  generic :: get => get_integer, get_logical, get_real32, get_real64, get_string
END TYPE Metadata_type

interface new_Metadata
  module procedure new_Metadata
end interface
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
TYPE :: HaloExchange_type

! Purpose :
! -------
!   *HaloExchange* : 

! Methods :
! -------
!   setup : Setup using arrays detailing proc and glb_idx, bounds and parbound
!   execute : Do the halo exchange

! Author :
! ------
!   17-Dec-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure :: setup => HaloExchange__setup
  procedure, private :: HaloExchange__execute_int32_r1
  procedure, private :: HaloExchange__execute_int32_r2
  procedure, private :: HaloExchange__execute_int32_r3
  procedure, private :: HaloExchange__execute_real32_r1
  procedure, private :: HaloExchange__execute_real32_r2
  procedure, private :: HaloExchange__execute_real32_r3
  procedure, private :: HaloExchange__execute_real64_r1
  procedure, private :: HaloExchange__execute_real64_r2
  procedure, private :: HaloExchange__execute_real64_r3
  generic :: execute => &
      & HaloExchange__execute_int32_r1, &
      & HaloExchange__execute_int32_r2, &
      & HaloExchange__execute_int32_r3, &
      & HaloExchange__execute_real32_r1, &
      & HaloExchange__execute_real32_r2, &
      & HaloExchange__execute_real32_r3, &
      & HaloExchange__execute_real64_r1, &
      & HaloExchange__execute_real64_r2, &
      & HaloExchange__execute_real64_r3
END TYPE HaloExchange_type
!------------------------------------------------------------------------------

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

INTERFACE view1d
  module procedure view1d_int32_rank2
  module procedure view1d_int32_rank3
  module procedure view1d_real32_rank2
  module procedure view1d_real32_rank3
  module procedure view1d_real64_rank2
  module procedure view1d_real64_rank3
  module procedure view1d_real64_rank4
end interface view1d

! =============================================================================
CONTAINS
! =============================================================================

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


! -----------------------------------------------------------------------------
! Helper functions

function c_to_f_string_str(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  character(len=:), allocatable :: mold
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  allocate(character(len=nchars) :: mold)
  str = transfer(s(1:nchars), mold)
end function c_to_f_string_str

function c_to_f_string_cptr(cptr) result(str)
  use iso_c_binding
  type(c_ptr), intent(in) :: cptr
  character(len=:), allocatable :: str
  character, dimension(:), pointer  :: s
  call C_F_POINTER ( cptr , s, (/MAX_STR_LEN/) )
  str = c_to_f_string_str(s)
end function c_to_f_string_cptr

function c_str(f_str)
  use, intrinsic :: iso_c_binding
  character(len=*), intent(in) :: f_str
  character(len=len_trim(f_str)+1) :: c_str
  c_str = trim(f_str) // c_null_char
end function c_str

function c_loc_int32(x)
  use iso_c_binding
  integer, target :: x
  type(c_ptr) :: c_loc_int32
  c_loc_int32 = C_LOC(x)
end function

function c_loc_real32(x)
  use iso_c_binding
  real(c_float), target :: x
  type(c_ptr) :: c_loc_real32
  c_loc_real32 = C_LOC(x)
end function

function c_loc_real64(x)
  use iso_c_binding
  real(c_double), target :: x
  type(c_ptr) :: c_loc_real64
  c_loc_real64 = C_LOC(x)
end function

function view1d_int32_rank2(array) result( view )
  integer, intent(in), target :: array(:,:)
  type(c_ptr) :: array_c_ptr
  integer, pointer :: view(:)
  array_c_ptr = c_loc_int32(array(1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_int32_rank2

function view1d_int32_rank3(array) result( view )
  integer, intent(in), target :: array(:,:,:)
  type(c_ptr) :: array_c_ptr
  integer, pointer :: view(:)
  array_c_ptr = c_loc_int32(array(1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_int32_rank3

function view1d_real32_rank2(array) result( view )
  real(c_float), intent(in), target :: array(:,:)
  type(c_ptr) :: array_c_ptr
  real(c_float), pointer :: view(:)
#ifndef  __GFORTRAN__
  if( .not. is_contiguous(array) ) then
    write(0,*) "ERROR: array is not contiguous in view1d"
    write(0,*) 'call abort()'
  end if
#endif
  array_c_ptr = c_loc_real32(array(1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real32_rank2

function view1d_real32_rank3(array) result( view )
  real(c_float), intent(in), target :: array(:,:,:)
  type(c_ptr) :: array_c_ptr
  real(c_float), pointer :: view(:)
  array_c_ptr = c_loc_real32(array(1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real32_rank3

function view1d_real64_rank2(array) result( view )
  real(c_double), intent(in), target :: array(:,:)
  type(c_ptr) :: array_c_ptr
  real(c_double), pointer :: view(:)
  array_c_ptr = c_loc_real64(array(1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real64_rank2

function view1d_real64_rank3(array) result( view )
  real(c_double), intent(in), target :: array(:,:,:)
  type(c_ptr) :: array_c_ptr
  real(c_double), pointer :: view(:)
  array_c_ptr = c_loc_real64(array(1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real64_rank3

function view1d_real64_rank4(array) result( view )
  real(c_double), intent(in), target :: array(:,:,:,:)
  type(c_ptr) :: array_c_ptr
  real(c_double), pointer :: view(:)
  array_c_ptr = c_loc_real64(array(1,1,1,1))
  call C_F_POINTER ( array_c_ptr , view , (/size(array)/) )
end function view1d_real64_rank4

! -----------------------------------------------------------------------------
! Mesh routines

function new_Mesh() result(mesh)
  type(Mesh_type) :: mesh
  mesh%object = atlas__Mesh__new()
end function new_Mesh

subroutine Mesh__add_function_space(this, function_space)
  class(Mesh_type), intent(inout) :: this
  type(FunctionSpace_type), intent(in) :: function_space
  call atlas__Mesh__add_function_space(this%object,function_space%object)
end subroutine Mesh__add_function_space

function Mesh__function_space(this,name) result(function_space)
  class(Mesh_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(FunctionSpace_type) :: function_space
  function_space%object = atlas__Mesh__function_space(this%object, c_str(name) )
  if( .not. C_associated(function_space%object) ) write(0,*) 'call abort()'
end function Mesh__function_space

subroutine Mesh__delete(this)
  type(Mesh_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call atlas__Mesh__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine Mesh__delete

! -----------------------------------------------------------------------------
! FunctionSpace routines

function new_FunctionSpace(name,shape_func,nb_nodes) result(function_space)
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_nodes
  type(FunctionSpace_type) :: function_space
  function_space%object = atlas__FunctionSpace__new(c_str(name),c_str(shape_func), &
    & (/nb_nodes,FIELD_NB_VARS/), 2 )
end function new_FunctionSpace

function new_PrismaticFunctionSpace(name,shape_func,nb_levels,nb_nodes) result(function_space)
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_levels
  integer, intent(in) :: nb_nodes
  type(FunctionSpace_type) :: function_space
  function_space%object = atlas__FunctionSpace__new(c_str(name),c_str(shape_func), &
    & (/nb_nodes,nb_levels,FIELD_NB_VARS/), 3 )
end function new_PrismaticFunctionSpace

subroutine FunctionSpace__delete(this)
  type(FunctionSpace_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call atlas__FunctionSpace__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine FunctionSpace__delete

subroutine FunctionSpace__create_field(this,name,nvars,kind)
  class(FunctionSpace_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: nvars
  integer, intent(in), optional :: kind
  if (present(kind)) then
    if (kind == KIND_REAL64) then
      call atlas__FunctionSpace__create_field_double(this%object,c_str(name),nvars)
    else if (kind == KIND_REAL32) then
      call atlas__FunctionSpace__create_field_float(this%object,c_str(name),nvars)
    else if (kind == KIND_INT32) then
      call atlas__FunctionSpace__create_field_int(this%object,c_str(name),nvars)
    else
      write(0,*) "Unsupported kind"
      write(0,*) 'call abort()'
    endif
  else if (wp == c_double) then
    call atlas__FunctionSpace__create_field_double(this%object,c_str(name),nvars)
  else if (wp == c_float) then
    call atlas__FunctionSpace__create_field_float(this%object,c_str(name),nvars)
  end if
end subroutine FunctionSpace__create_field


subroutine FunctionSpace__remove_field(this,name)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  call atlas__FunctionSpace__remove_field(this%object,c_str(name))
end subroutine FunctionSpace__remove_field

function FunctionSpace__name(this) result(name)
  class(FunctionSpace_type), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = atlas__FunctionSpace__name(this%object)
  name = c_to_f_string_cptr(name_c_str)
end function FunctionSpace__name

function FunctionSpace__dof(this) result(dof)
  class(FunctionSpace_type), intent(in) :: this
  integer :: dof
  dof = atlas__FunctionSpace__dof(this%object)
end function FunctionSpace__dof

function FunctionSpace__glb_dof(this) result(glb_dof)
  class(FunctionSpace_type), intent(in) :: this
  integer :: glb_dof
  glb_dof = atlas__FunctionSpace__glb_dof(this%object)
end function FunctionSpace__glb_dof

function FunctionSpace__bounds(this) result(bounds)
  class(FunctionSpace_type), intent(in) :: this
  integer, pointer :: bounds(:)
  type(c_ptr) :: bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__FunctionSpace__boundsf(this%object, bounds_c_ptr, field_rank)
  call C_F_POINTER ( bounds_c_ptr , bounds , (/field_rank/) )
end function FunctionSpace__bounds

function FunctionSpace__field(this,name) result(field)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(Field_type) :: field
  field%object = atlas__FunctionSpace__field(this%object, c_str(name) )
  if( .not. C_associated(field%object) ) write(0,*) 'call abort()'
end function FunctionSpace__field

function FunctionSpace__has_field(this,name) result(flag)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FunctionSpace__has_field(this%object, c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FunctionSpace__has_field

subroutine FunctionSpace__parallelise(this)
  class(FunctionSpace_type), intent(in) :: this
  call atlas__FunctionSpace__parallelise(this%object)
end subroutine FunctionSpace__parallelise

function FunctionSpace__get_halo_exchange(this) result(halo_exchange)
  class(FunctionSpace_type), intent(in) :: this
  type(HaloExchange_type) :: halo_exchange
  halo_exchange%object = atlas__FunctionSpace__halo_exchange( this%object )
end function FunctionSpace__get_halo_exchange

subroutine FunctionSpace__halo_exchange_int32_r1(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_int( this%object, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r1
subroutine FunctionSpace__halo_exchange_int32_r2(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r2
subroutine FunctionSpace__halo_exchange_int32_r3(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_int( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_int32_r3

subroutine FunctionSpace__halo_exchange_real32_r1(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_float( this%object, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r1
subroutine FunctionSpace__halo_exchange_real32_r2(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r2
subroutine FunctionSpace__halo_exchange_real32_r3(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_float( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real32_r3

subroutine FunctionSpace__halo_exchange_real64_r1(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__halo_exchange_double( this%object, field_data, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r1
subroutine FunctionSpace__halo_exchange_real64_r2(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r2
subroutine FunctionSpace__halo_exchange_real64_r3(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r3
subroutine FunctionSpace__halo_exchange_real64_r4(this, field_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:,:)
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__FunctionSpace__halo_exchange_double( this%object, view, size(field_data) )
end subroutine FunctionSpace__halo_exchange_real64_r4


subroutine FunctionSpace__gather_real32_r1(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(in) :: field_data(:)
  real(c_float), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_float( this%object, field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r1


subroutine FunctionSpace__gather_real32_r2(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_float( this%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r2


subroutine FunctionSpace__gather_real32_r3(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_float), intent(in) :: field_data(:,:,:)
  real(c_float), intent(inout) :: glbfield_data(:,:,:)
  real(c_float), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_float( this%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real32_r3

subroutine FunctionSpace__gather_real64_r1(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(in) :: field_data(:)
  real(c_double), intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_double( this%object, field_data, size(field_data), &
                                          & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r1


subroutine FunctionSpace__gather_real64_r2(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_double( this%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r2


subroutine FunctionSpace__gather_real64_r3(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  real(c_double), intent(in) :: field_data(:,:,:)
  real(c_double), intent(inout) :: glbfield_data(:,:,:)
  real(c_double), pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_double( this%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_real64_r3

subroutine FunctionSpace__gather_int32_r1(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(in) :: field_data(:)
  integer, intent(inout) :: glbfield_data(:)
#ifndef  __GFORTRAN__
  if (.not. is_contiguous(field_data) ) then
    write(0,*) "ERROR: field_data is not contiguous"
    write(0,*) 'call abort()'
  end if
#endif
  call atlas__FunctionSpace__gather_int( this%object, field_data, size(field_data), &
                                       & glbfield_data, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r1

subroutine FunctionSpace__gather_int32_r2(this, field_data, glbfield_data)
  class(FunctionSpace_type), intent(in) :: this
  integer, intent(in) :: field_data(:,:)
  integer, intent(inout) :: glbfield_data(:,:)
  integer, pointer :: view(:), glbview(:)
  view => view1d(field_data)
  glbview => view
  if( size(glbfield_data) /= 0 ) then
    glbview => view1d(glbfield_data)
  end if
  call atlas__FunctionSpace__gather_int( this%object, view, size(field_data), &
                                          & glbview, size(glbfield_data) )
end subroutine FunctionSpace__gather_int32_r2




! ------------------------------------------------------------------------------
! HaloExchange routines

function new_HaloExchange() result(halo_exchange)
  type(HaloExchange_type) :: halo_exchange
  halo_exchange%object = atlas__HaloExchange__new()
end function new_HaloExchange

subroutine HaloExchange__delete(this)
  type(HaloExchange_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call atlas__HaloExchange__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine HaloExchange__delete

subroutine HaloExchange__setup(this, part, remote_idx)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(in) :: part(:)
  integer, intent(in) :: remote_idx(:)
  call atlas__HaloExchange__setup( this%object, part, remote_idx, size(part) )
end subroutine HaloExchange__setup


subroutine HaloExchange__execute_int32_r1(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__HaloExchange__execute_int( this%object, field_data, 1 )
  else
    call atlas__HaloExchange__execute_int( this%object, field_data, nb_vars )
  end if
end subroutine HaloExchange__execute_int32_r1


subroutine HaloExchange__execute_int32_r2(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:)
  integer, intent(in) :: nb_vars
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_int( this%object, view, nb_vars )
end subroutine HaloExchange__execute_int32_r2


subroutine HaloExchange__execute_int32_r3(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  integer, intent(inout) :: field_data(:,:,:)
  integer, intent(in) :: nb_vars
  integer, pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_int( this%object, view, nb_vars )
end subroutine HaloExchange__execute_int32_r3

subroutine HaloExchange__execute_real32_r1(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__HaloExchange__execute_float( this%object, field_data, 1 )
  else
    call atlas__HaloExchange__execute_float( this%object, field_data, nb_vars )
  end if
end subroutine HaloExchange__execute_real32_r1
subroutine HaloExchange__execute_real32_r2(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:)
  integer, intent(in) :: nb_vars
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_float( this%object, view, nb_vars )
end subroutine HaloExchange__execute_real32_r2
subroutine HaloExchange__execute_real32_r3(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_float), intent(inout) :: field_data(:,:,:)
  integer, intent(in) :: nb_vars
  real(c_float), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_float( this%object, view, nb_vars )
end subroutine HaloExchange__execute_real32_r3

subroutine HaloExchange__execute_real64_r1(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:)
  integer, optional, intent(in) :: nb_vars
  if (.not. present(nb_vars) ) then
    call atlas__HaloExchange__execute_double( this%object, field_data, 1 )
  else
    call atlas__HaloExchange__execute_double( this%object, field_data, nb_vars )
  end if
end subroutine HaloExchange__execute_real64_r1
subroutine HaloExchange__execute_real64_r2(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:)
  integer, intent(in) :: nb_vars
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_double( this%object, view, nb_vars )
end subroutine HaloExchange__execute_real64_r2
subroutine HaloExchange__execute_real64_r3(this, field_data, nb_vars)
  class(HaloExchange_type), intent(in) :: this
  real(c_double), intent(inout) :: field_data(:,:,:)
  integer, intent(in) :: nb_vars
  real(c_double), pointer :: view(:)
  view => view1d(field_data)
  call atlas__HaloExchange__execute_double( this%object, view, nb_vars )
end subroutine HaloExchange__execute_real64_r3



! -----------------------------------------------------------------------------
! Field routines

function Field__name(this) result(field_name)
  class(Field_type), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = atlas__Field__name(this%object)
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

function Field__data_type(this) result(field_data_type)
  class(Field_type), intent(in) :: this
  character(len=:), allocatable :: field_data_type
  type(c_ptr) :: field_data_type_c_str
  field_data_type_c_str = atlas__Field__data_type(this%object)
  field_data_type = c_to_f_string_cptr(field_data_type_c_str)
end function Field__data_type

function Field__nb_vars(this) result(nb_vars)
  class(Field_type), intent(in) :: this
  integer :: nb_vars
  nb_vars = atlas__Field__nb_vars(this%object)
end function Field__nb_vars

function Field__metadata(this) result(metadata)
  class(Field_type), intent(in) :: this
  type(metadata_type) :: Metadata
  metadata%object = atlas__Field__metadata(this%object)
end function Field__metadata

function Field__function_space(this) result(function_space)
  class(Field_type), intent(in) :: this
  type(FunctionSpace_type) :: function_space
  function_space%object = atlas__Field__function_space(this%object)
end function Field__function_space

subroutine Field__access_data1_integer(this, field) 
  class(Field_type), intent(in) :: this
  integer, pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_boundsf_int(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_integer

subroutine Field__access_data2_integer(this, field) 
  class(Field_type), intent(in) :: this
  integer, pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_int(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_integer

subroutine Field__access_data3_integer(this, field) 
  class(Field_type), intent(in) :: this
  integer, pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_int(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_integer

subroutine Field__access_data1_real32(this, field) 
  class(Field_type), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_boundsf_float(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_real32

subroutine Field__access_data2_real32(this, field) 
  class(Field_type), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_float(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_real32

subroutine Field__access_data3_real32(this, field) 
  class(Field_type), intent(in) :: this
  real(c_float), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_float(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_real32

subroutine Field__access_data1_real64(this, field) 
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call atlas__Field__data_boundsf_double(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end subroutine Field__access_data1_real64

subroutine Field__access_data2_real64(this, field) 
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_double(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
end subroutine Field__access_data2_real64

subroutine Field__access_data3_real64(this, field) 
  class(Field_type), intent(in) :: this
  real(c_double), pointer, intent(out) :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call atlas__Field__data_boundsf_double(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-2:field_rank) )
end subroutine Field__access_data3_real64



function Field__data1_wp(this) result(field)
  class(Field_type), intent(in) :: this
  real(wp), pointer :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  if( wp == c_double ) then
    call atlas__Field__data_boundsf_double(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  else if (wp == c_float ) then
    call atlas__Field__data_boundsf_float(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  end if
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end function Field__data1_wp

function Field__data2_wp(this) result(field)
  class(Field_type), intent(in) :: this
  real(wp), pointer :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  if( wp == c_double ) then
    call atlas__Field__data_boundsf_double(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  else if (wp == c_float ) then
    call atlas__Field__data_boundsf_float(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  end if
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  if (size(field_bounds) < 2) then
    write(0,*) "Cannot access field """,this%name(),""" with rank",field_rank," as rank 2"
    write(0,*) 'call abort()'
  end if
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , field_bounds(field_rank-1:field_rank) )
  if( size(field) /= field_size ) then
    write(0,*) "Requested bounds of field ", this%name(), "[", field_bounds(1:2), &
     & "] do not cover the entire field of size ", field_size
    write(0,*) 'call abort()'
  end if
end function Field__data2_wp

function Field__data3_wp(this) result(field)
  class(Field_type), intent(in) :: this
  real(wp), pointer :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  if( wp == c_double ) then
    call atlas__Field__data_boundsf_double(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  else if (wp == c_float ) then
    call atlas__Field__data_boundsf_float(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  end if
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  if (size(field_bounds) < 3) then
    write(0,*) "Cannot access field """,this%name(),""" with rank",field_rank," as rank 3"
    write(0,*) 'call abort()'
  end if
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , field_bounds(1:3) )
  if( size(field) /= field_size ) then
    write(0,*) "Requested bounds of field ", field_bounds(1:3), " do not cover the entire field of size ", field_size
    write(0,*) 'call abort()'
  end if
end function Field__data3_wp

! -----------------------------------------------------------------------------
! FieldSet routines

function new_FieldSet(name) result(fieldset)
  character(len=*), intent(in) :: name
  type(FieldSet_type) :: fieldset
  fieldset%object = atlas__FieldSet__new( c_str(name) )
end function new_FieldSet

subroutine FieldSet__delete(this)
  type(FieldSet_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call atlas__FieldSet__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine FieldSet__delete

subroutine FieldSet__add_field(this,field)
  class(FieldSet_type), intent(in) :: this
  type(Field_type), intent(in) :: field
  call atlas__FieldSet__add_field(this%object, field%object)
end subroutine FieldSet__add_field

function FieldSet__has_field(this,name) result(flag)
  class(FieldSet_type), intent(in) :: this
  character(len=*), intent(in) :: name
  logical :: flag
  integer :: rc
  rc = atlas__FieldSet__has_field(this%object, c_str(name))
  if( rc == 0 ) then
    flag = .False.
  else
    flag = .True.
  end if
end function FieldSet__has_field

function FieldSet__size(this) result(nb_fields)
  class(FieldSet_type), intent(in) :: this
  integer :: nb_fields
  nb_fields = atlas__FieldSet__size(this%object)
end function FieldSet__size

function FieldSet__field_by_name(this,name) result(field)
  class(FieldSet_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(Field_type) :: field
  field%object = atlas__FieldSet__field_by_name(this%object, c_str(name) )
end function FieldSet__field_by_name

function FieldSet__field_by_idx(this,idx) result(field)
  class(FieldSet_type), intent(in) :: this
  integer, intent(in) :: idx
  type(Field_type) :: field
  field%object = atlas__FieldSet__field_by_idx(this%object, idx-1) ! C index
end function FieldSet__field_by_idx

subroutine FieldSet__fields(this,fields)
  class(FieldSet_type), intent(in) :: this
  type(Field_type), allocatable, intent(out) :: fields(:)

  type(c_ptr), pointer :: fields_ptr(:)
  type(c_ptr) :: fields_cptr
  integer :: nb_fields, jfield
  call atlas__FieldSet__fields(this%object, fields_cptr, nb_fields)
  call c_f_pointer( fields_cptr, fields_ptr, (/nb_fields/) )
  allocate( fields(nb_fields) )
  do jfield=1,nb_fields
    fields(jfield)%object = fields_ptr(jfield)
  end do
end subroutine FieldSet__fields


! -----------------------------------------------------------------------------
! Metadata routines

function new_Metadata() result(metadata)
  type(Metadata_type) :: metadata
  metadata%object = atlas__Metadata__new()
end function new_Metadata

subroutine Metadata__delete(this)
  type(Metadata_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call atlas__Metadata__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine Metadata__delete

subroutine Metadata__add_logical(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  logical, intent(in) :: value
  integer :: value_int
  if( value ) then 
    value_int = 1
  else
    value_int = 0
  end if
  call atlas__Metadata__add_int(this%object, c_str(name), value_int )
end subroutine Metadata__add_logical

subroutine Metadata__add_integer(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call atlas__Metadata__add_int(this%object, c_str(name), value)
end subroutine Metadata__add_integer

subroutine Metadata__add_real32(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(in) :: value
  call atlas__Metadata__add_float(this%object, c_str(name) ,value)
end subroutine Metadata__add_real32

subroutine Metadata__add_real64(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(in) :: value
  call atlas__Metadata__add_double(this%object, c_str(name) ,value)
end subroutine Metadata__add_real64

subroutine Metadata__add_string(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call atlas__Metadata__add_string(this%object, c_str(name) , c_str(value) )
end subroutine Metadata__add_string

subroutine Metadata__get_logical(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(out) :: value
  integer :: value_int
  value_int = atlas__Metadata__get_int(this%object,c_str(name) )
  if (value_int > 0) then
    value = .True.
  else 
    value = .False.
  end if
end subroutine Metadata__get_logical

subroutine Metadata__get_integer(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  integer, intent(out) :: value
  value = atlas__Metadata__get_int(this%object, c_str(name) )
end subroutine Metadata__get_integer

subroutine Metadata__get_real32(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_float), intent(out) :: value
  value = atlas__Metadata__get_float(this%object, c_str(name) )
end subroutine Metadata__get_real32

subroutine Metadata__get_real64(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  real(c_double), intent(out) :: value
  value = atlas__Metadata__get_double(this%object, c_str(name) )
end subroutine Metadata__get_real64

subroutine Metadata__get_string(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(out) :: value
  type(c_ptr) :: value_c_str
  value_c_str = atlas__Metadata__get_string(this%object, c_str(name) )
  value = c_to_f_string_cptr(value_c_str)
end subroutine Metadata__get_string

! -----------------------------------------------------------------------------

function atlas_read_gmsh(filename) result(mesh)
  character(len=*), intent(in) :: filename
  type(Mesh_type) :: mesh
  mesh%object = atlas__read_gmsh(c_str(filename))
end function atlas_read_gmsh

subroutine atlas_write_gmsh(mesh,filename)
  type(Mesh_type), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_gmsh(mesh%object,c_str(filename))
end subroutine atlas_write_gmsh

subroutine atlas_build_parallel_fields(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_parallel_fields(mesh%object)
end subroutine atlas_build_parallel_fields

subroutine atlas_build_periodic_boundaries(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_periodic_boundaries(mesh%object)
end subroutine atlas_build_periodic_boundaries

subroutine atlas_build_halo(mesh,nelems)
  type(Mesh_type), intent(inout) :: mesh
  integer, intent(in) :: nelems
  call atlas__build_halo(mesh%object,nelems)
end subroutine atlas_build_halo

subroutine atlas_build_edges(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_edges(mesh%object)
end subroutine atlas_build_edges

subroutine atlas_build_pole_edges(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_pole_edges(mesh%object)
end subroutine atlas_build_pole_edges

subroutine atlas_build_dual_mesh(mesh)
  type(Mesh_type), intent(inout) :: mesh
  call atlas__build_dual_mesh(mesh%object)
end subroutine atlas_build_dual_mesh

! -----------------------------------------------------------------------------

end module atlas_module
