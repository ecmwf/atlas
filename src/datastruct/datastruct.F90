
module datastruct

! Purpose :
! -------
!    *datastruct* : Low-level API to manipulate meshes in a
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

! Interfaces :
! ----------

! Author :
! ------
!   20-Nov-2013 Willem Deconinck    *ECMWF*

!------------------------------------------------------------------------------
use, intrinsic :: iso_c_binding
use field_c_binding
use fieldset_c_binding
use functionspace_c_binding
use mesh_c_binding
use metadata_c_binding
implicit none

integer, private, parameter :: MAX_STR_LEN = 255

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
!   create_field : Create a new field in this function space with given name
!   remove_field : Remove a field with given name
!   field : Access to a field with given name

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
private
  type(c_ptr) :: object = C_NULL_ptr
contains
  procedure :: name => FunctionSpace__name
  procedure :: create_field => FunctionSpace__create_field
  procedure :: remove_field => FunctionSpace__remove_field
  procedure :: field => FunctionSpace__field
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
  procedure :: Metadata => Field__Metadata
  procedure :: data1 => Field__data1
  procedure :: data2 => Field__data2
  procedure :: data3 => Field__data3
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
  procedure, private :: add_real4 => Metadata__add_real4
  procedure, private :: add_real8 => Metadata__add_real8
  procedure, private :: add_string => Metadata__add_string
  generic :: add => add_logical, add_integer, add_real4, add_real8, add_string
  generic :: set => add_logical, add_integer, add_real4, add_real8, add_string
  procedure :: get_integer => Metadata__get_integer
  procedure :: get_logical => Metadata__get_logical
  procedure :: get_real4 => Metadata__get_real4
  procedure :: get_real8 => Metadata__get_real8
  procedure :: get_string => Metadata__get_string
  generic :: get => get_integer, get_logical, get_real4, get_real8, get_string
END TYPE Metadata_type
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
end interface delete

! =============================================================================
CONTAINS
! =============================================================================

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
  use iso_c_binding
  character(len=*), intent(in) :: f_str
  character(len=len_trim(f_str)+1) :: c_str
  c_str = trim(f_str) // c_null_char
end function c_str

! -----------------------------------------------------------------------------
! Mesh routines

function new_Mesh() result(mesh)
  type(Mesh_type) :: mesh
  mesh%object = ecmwf__Mesh__new()
end function new_Mesh

subroutine Mesh__add_function_space(this, function_space)
  class(Mesh_type), intent(inout) :: this
  type(FunctionSpace_type), intent(in) :: function_space
  call ecmwf__Mesh__add_function_space(this%object,function_space%object)
end subroutine Mesh__add_function_space

function Mesh__function_space(this,name) result(function_space)
  class(Mesh_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(FunctionSpace_type) :: function_space
  function_space%object = ecmwf__Mesh__function_space(this%object, c_str(name) )
  if( .not. C_associated(function_space%object) ) call abort()
end function Mesh__function_space

subroutine Mesh__delete(this)
  type(Mesh_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call ecmwf__Mesh__delete(this%object)
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
  function_space%object = ecmwf__FunctionSpace__new(c_str(name),c_str(shape_func),nb_nodes)
end function new_FunctionSpace

function new_PrismaticFunctionSpace(name,shape_func,nb_levels,nb_nodes) result(function_space)
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: shape_func
  integer, intent(in) :: nb_levels
  integer, intent(in) :: nb_nodes
  type(FunctionSpace_type) :: function_space
  function_space%object = ecmwf__PrismaticFunctionSpace__new(c_str(name),c_str(shape_func),nb_levels,nb_nodes)
end function new_PrismaticFunctionSpace

subroutine FunctionSpace__delete(this)
  type(FunctionSpace_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call ecmwf__FunctionSpace__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine FunctionSpace__delete

subroutine FunctionSpace__create_field(this,name,size)
  class(FunctionSpace_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: size
  call ecmwf__FunctionSpace__create_field(this%object,c_str(name),size)
end subroutine FunctionSpace__create_field

subroutine FunctionSpace__remove_field(this,name)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  call ecmwf__FunctionSpace__remove_field(this%object,c_str(name))
end subroutine FunctionSpace__remove_field

function FunctionSpace__name(this) result(name)
  class(FunctionSpace_type), intent(in) :: this
  character(len=:), allocatable :: name
  type(c_ptr) :: name_c_str
  name_c_str = ecmwf__FunctionSpace__name(this%object)
  name = c_to_f_string_cptr(name_c_str)
end function FunctionSpace__name


function FunctionSpace__field(this,name) result(field)
  class(FunctionSpace_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(Field_type) :: field
  field%object = ecmwf__FunctionSpace__field(this%object, c_str(name) )
  if( .not. C_associated(field%object) ) call abort()
end function FunctionSpace__field

! -----------------------------------------------------------------------------
! Field routines

function Field__name(this) result(field_name)
  class(Field_type), intent(in) :: this
  character(len=:), allocatable :: field_name
  type(c_ptr) :: field_name_c_str
  field_name_c_str = ecmwf__Field__name(this%object)
  field_name = c_to_f_string_cptr(field_name_c_str)
end function Field__name

function Field__Metadata(this) result(Metadata)
  class(Field_type), intent(in) :: this
  type(Metadata_type) :: Metadata
  Metadata%object = ecmwf__Field__Metadata(this%object)
end function Field__Metadata

function Field__data1(this) result(field)
  class(Field_type), intent(in) :: this
  real(kind=8), pointer :: field(:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  integer :: field_size, jbound
  call ecmwf__Field__data(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  field_size = 1
  do jbound=1,field_rank
    field_size = field_size * field_bounds(jbound)
  end do
  call C_F_POINTER ( field_c_ptr , field , (/field_size/) )
end function Field__data1

function Field__data2(this) result(field)
  class(Field_type), intent(in) :: this
  real(kind=8), pointer :: field(:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call ecmwf__Field__data(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds(1:2) )
end function Field__data2

function Field__data3(this) result(field)
  class(Field_type), intent(in) :: this
  real(kind=8), pointer :: field(:,:,:)
  integer, pointer :: field_bounds(:)
  type(c_ptr) :: field_c_ptr
  type(c_ptr) :: field_bounds_c_ptr
  integer(c_int) :: field_rank
  call ecmwf__Field__data(this%object, field_c_ptr, field_bounds_c_ptr, field_rank)
  call C_F_POINTER ( field_bounds_c_ptr , field_bounds , (/field_rank/) )
  call C_F_POINTER ( field_c_ptr , field , field_bounds )
end function Field__data3

! -----------------------------------------------------------------------------
! FieldSet routines

function new_FieldSet(name) result(fieldset)
  character(len=*), intent(in) :: name
  type(FieldSet_type) :: fieldset
  fieldset%object = ecmwf__FieldSet__new( c_str(name) )
end function new_FieldSet

subroutine FieldSet__delete(this)
  type(FieldSet_type), intent(inout) :: this
  if ( c_associated(this%object) ) then
    call ecmwf__FieldSet__delete(this%object)
  end if
  this%object = C_NULL_ptr
end subroutine FieldSet__delete

subroutine FieldSet__add_field(this,field)
  class(FieldSet_type), intent(in) :: this
  type(Field_type), intent(in) :: field
  call ecmwf__FieldSet__add_field(this%object, field%object)
end subroutine FieldSet__add_field

function FieldSet__size(this) result(nb_fields)
  class(FieldSet_type), intent(in) :: this
  integer :: nb_fields
  nb_fields = ecmwf__FieldSet__size(this%object)
end function FieldSet__size

function FieldSet__field_by_name(this,name) result(field)
  class(FieldSet_type), intent(in) :: this
  character(len=*), intent(in) :: name
  type(Field_type) :: field
  field%object = ecmwf__FieldSet__field_by_name(this%object, c_str(name) )
end function FieldSet__field_by_name

function FieldSet__field_by_idx(this,idx) result(field)
  class(FieldSet_type), intent(in) :: this
  integer, intent(in) :: idx
  type(Field_type) :: field
  field%object = ecmwf__FieldSet__field_by_idx(this%object, idx-1) ! C index
end function FieldSet__field_by_idx

subroutine FieldSet__fields(this,fields)
  class(FieldSet_type), intent(in) :: this
  type(Field_type), allocatable, intent(out) :: fields(:)

  type(c_ptr), pointer :: fields_ptr(:)
  type(c_ptr) :: fields_cptr
  integer :: nb_fields, jfield
  call ecmwf__FieldSet__fields(this%object, fields_cptr, nb_fields)
  call c_f_pointer( fields_cptr, fields_ptr, (/nb_fields/) )
  allocate( fields(nb_fields) )
  do jfield=1,nb_fields
    fields(jfield)%object = fields_ptr(jfield)
  end do
end subroutine FieldSet__fields

! -----------------------------------------------------------------------------
! Metadata routines

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
  call ecmwf__Metadata__add_int(this%object, c_str(name), value_int )
end subroutine Metadata__add_logical

subroutine Metadata__add_integer(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  integer, intent(in) :: value
  call ecmwf__Metadata__add_int(this%object, c_str(name), value)
end subroutine Metadata__add_integer

subroutine Metadata__add_real4(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(kind=4), intent(in) :: value
  call ecmwf__Metadata__add_float(this%object, c_str(name) ,value)
end subroutine Metadata__add_real4

subroutine Metadata__add_real8(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  real(kind=8), intent(in) :: value
  call ecmwf__Metadata__add_double(this%object, c_str(name) ,value)
end subroutine Metadata__add_real8

subroutine Metadata__add_string(this, name, value)
  class(Metadata_type), intent(inout) :: this
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: value
  call ecmwf__Metadata__add_string(this%object, c_str(name) , c_str(value) )
end subroutine Metadata__add_string

subroutine Metadata__get_logical(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  logical, intent(out) :: value
  integer :: value_int
  value_int = ecmwf__Metadata__get_int(this%object,c_str(name) )
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
  value = ecmwf__Metadata__get_int(this%object, c_str(name) )
end subroutine Metadata__get_integer

subroutine Metadata__get_real4(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  real(kind=4), intent(out) :: value
  value = ecmwf__Metadata__get_float(this%object, c_str(name) )
end subroutine Metadata__get_real4

subroutine Metadata__get_real8(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  real(kind=8), intent(out) :: value
  value = ecmwf__Metadata__get_double(this%object, c_str(name) )
end subroutine Metadata__get_real8

subroutine Metadata__get_string(this, name, value)
  class(Metadata_type), intent(in) :: this
  character(len=*), intent(in) :: name
  character(len=:), allocatable, intent(out) :: value
  type(c_ptr) :: value_c_str
  value_c_str = ecmwf__Metadata__get_string(this%object, c_str(name) )
  value = c_to_f_string_cptr(value_c_str)
end subroutine Metadata__get_string

! -----------------------------------------------------------------------------

end module datastruct