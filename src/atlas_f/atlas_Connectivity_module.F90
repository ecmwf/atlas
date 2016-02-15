
module atlas_connectivity_module

use iso_c_binding, only : c_funptr, c_ptr, c_loc, c_f_pointer, c_f_procpointer, c_funloc, c_int, c_size_t
use atlas_refcounted_module, only : atlas_refcounted
implicit none

private :: c_funptr
private :: c_ptr
private :: c_loc
private :: c_f_pointer
private :: c_funloc
private :: c_int
private :: c_size_t
private :: atlas_refcounted

public :: atlas_Connectivity
public :: atlas_MultiBlockConnectivity
public :: atlas_BlockConnectivity


private

!-----------------------------
! atlas_Connectivity         !
!-----------------------------

type, extends(atlas_refcounted) :: atlas_Connectivity

! Public members
  type( atlas_ConnectivityAccess ), public, pointer :: access => null()

contains

! Public methods
  procedure, private :: value_args_int       => atlas_Connectivity__value_args_int
  procedure, private :: value_args_size_t    => atlas_Connectivity__value_args_size_t
  generic, public :: value => value_args_int, value_args_size_t
  procedure, public :: rows     => atlas_Connectivity__rows
  procedure, public :: cols     => atlas_Connectivity__cols
  procedure, public :: maxcols  => atlas_Connectivity__maxcols
  procedure, public :: mincols  => atlas_Connectivity__mincols
  procedure, public :: padded_data     => atlas_Connectivity__padded_data
  procedure, public :: data     => atlas_Connectivity__data
  procedure, public :: row      => atlas_Connectivity__row
  procedure, public :: copy     => atlas_Connectivity__copy
  procedure, public :: delete   => atlas_Connectivity__delete
  procedure, public :: missing_value  => atlas_Connectivity__missing_value
  procedure, private :: add_values_args_int => atlas_Connectivity__add_values_args_int
  procedure, private :: add_values_args_size_t => atlas_Connectivity__add_values_args_size_t
  procedure, private :: add_missing_args_int => atlas_Connectivity__add_missing_args_int
  procedure, private :: add_missing_args_size_t => atlas_Connectivity__add_missing_args_size_t
  generic, public :: add => add_values_args_int, add_values_args_size_t, add_missing_args_int, add_missing_args_size_t
end type

!---------------------------------------
! atlas_MultiBlockConnectivity         !
!---------------------------------------

type, extends(atlas_Connectivity) :: atlas_MultiBlockConnectivity
contains
! Public methods
  procedure, public :: blocks   => atlas_MultiBlockConnectivity__blocks
  procedure, public :: block    => atlas_MultiBlockConnectivity__block
end type

!----------------------------!
! atlas_BlockConnectivity    !
!----------------------------!

type, extends(atlas_refcounted) :: atlas_BlockConnectivity
contains
  procedure, public :: copy     => atlas_BlockConnectivity__copy
  procedure, public :: delete   => atlas_BlockConnectivity__delete
  procedure, public :: rows     => atlas_BlockConnectivity__rows
  procedure, public :: cols     => atlas_BlockConnectivity__cols
  procedure, public :: data     => atlas_BlockConnectivity__data
  procedure, public :: missing_value  => atlas_BlockConnectivity__missing_value
end type

interface atlas_Connectivity
  module procedure Connectivity_cptr
  module procedure Connectivity_constructor
end interface

interface atlas_MultiBlockConnectivity
  module procedure MultiBlockConnectivity_cptr
  module procedure MultiBlockConnectivity_constructor
end interface

interface atlas_BlockConnectivity
  module procedure BlockConnectivity_cptr
end interface


!-------------------------------
! Helper types                 !
!-------------------------------

type :: atlas_ConnectivityAccessRow
  integer, pointer :: col(:)
contains
end type

type :: atlas_ConnectivityAccess
  integer(c_int),    private, pointer :: values_(:) => null()
  integer(c_size_t), private, pointer :: displs_(:) => null()
  integer(c_size_t), public,  pointer :: cols(:)    => null()
  type(atlas_ConnectivityAccessRow), public, pointer :: row(:) => null()
  integer(c_size_t), private :: rows_
  integer, private, pointer :: padded_(:,:) => null()
  integer(c_size_t), private :: maxcols_, mincols_
  integer(c_int), private :: missing_value_
  type(c_ptr), private :: connectivity_ptr
contains
  procedure, public :: rows  => access_rows
  procedure, public :: value => access_value
end type

! ----------------------------------------------------------------------------
interface
    subroutine atlas__connectivity__register_update(This,funptr,ptr) bind(c,name="atlas__connectivity__register_update")
        import :: c_funptr
        import :: c_ptr
        type(c_ptr),    value :: This
        type(c_funptr), value :: funptr
        type(c_ptr),    value :: ptr
    end subroutine

    subroutine atlas__connectivity__register_delete(This,funptr,ptr) bind(c,name="atlas__connectivity__register_delete")
        import :: c_funptr
        import :: c_ptr
        type(c_ptr),    value :: This
        type(c_funptr), value :: funptr
        type(c_ptr),    value :: ptr
    end subroutine

    function atlas__connectivity__ctxt_update(This,ptr) bind(c,name="atlas__connectivity__ctxt_update")
        import ::c_ptr
        import :: c_int
        integer(c_int) :: atlas__connectivity__ctxt_update
        type(c_ptr),    value :: This
        type(c_ptr) :: ptr
    end function

    function atlas__connectivity__ctxt_delete(This,ptr) bind(c,name="atlas__connectivity__ctxt_delete")
        import ::c_ptr
        import :: c_int
        integer(c_int) :: atlas__connectivity__ctxt_delete
        type(c_ptr),    value :: This
        type(c_ptr) :: ptr
    end function

end interface

contains

function Connectivity_cptr(cptr) result(this)
  use atlas_connectivity_c_binding
  type(atlas_Connectivity) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
  call setup_access(this)
end function

function Connectivity_constructor() result(this)
  use atlas_connectivity_c_binding
  type(atlas_Connectivity) :: this
  this = Connectivity_cptr( atlas__Connectivity__create() )
  call this%return()
end function

subroutine atlas_Connectivity__delete(this)
  use atlas_connectivity_c_binding
  class(atlas_Connectivity), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Connectivity__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_Connectivity__copy(this,obj_in)
  class(atlas_Connectivity), intent(inout) :: this
  class(atlas_RefCounted),   target, intent(in) :: obj_in
  class(atlas_Connectivity), pointer :: obj_in_cast
  select type (ptr => obj_in)
    type is (atlas_Connectivity)
      obj_in_cast => ptr
      if( associated( obj_in_cast%access ) ) this%access => obj_in_cast%access
    type is (atlas_MultiBlockConnectivity)
      obj_in_cast => ptr
      if( associated( obj_in_cast%access ) ) this%access => obj_in_cast%access
  end select
end subroutine

pure function access_value(this,c,r) result(val)
  integer(c_int) :: val
  class(atlas_ConnectivityAccess), intent(in) :: this
  integer(c_size_t), intent(in) :: r,c
  val = this%values_(c+this%displs_(r))
end function

pure function access_rows(this) result(val)
  integer(c_int) :: val
  class(atlas_ConnectivityAccess), intent(in) :: this
  val = this%rows_
end function

pure function atlas_Connectivity__value_args_size_t(this,c,r) result(val)
  integer(c_int) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_size_t), intent(in) :: r,c
  val = this%access%values_(c+this%access%displs_(r))
end function

pure function atlas_Connectivity__value_args_int(this,c,r) result(val)
  integer(c_int) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: r,c
  val = this%access%values_(c+this%access%displs_(r))
end function

pure function atlas_Connectivity__rows(this) result(val)
  integer(c_size_t) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%rows_
end function

function atlas_Connectivity__missing_value(this) result(val)
  use atlas_connectivity_c_binding
  integer(c_int) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = atlas__Connectivity__missing_value(this%c_ptr())
end function

pure function atlas_Connectivity__cols(this,r) result(val)
  integer(c_size_t) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_size_t), intent(in) :: r
 val = this%access%cols(r)
end function

pure function atlas_Connectivity__mincols(this) result(val)
  integer(c_size_t) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%mincols_
end function

pure function atlas_Connectivity__maxcols(this) result(val)
  integer(c_size_t) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%maxcols_
end function

subroutine atlas_Connectivity__padded_data(this, padded, cols)
  class(atlas_Connectivity), intent(inout) :: this
  integer(c_int), pointer, intent(inout) :: padded(:,:)
  integer(c_size_t), pointer, intent(out), optional :: cols(:)
  if( .not. associated(this%access%padded_) ) call update_padded(this%access)
  padded => this%access%padded_
  if( present(cols) ) cols => this%access%cols
end subroutine

subroutine atlas_Connectivity__data(this, data, ncols)
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), pointer, intent(out) :: data(:,:)
  integer(c_int), intent(out), optional :: ncols
  integer(c_int) :: maxcols
  maxcols = this%maxcols()
  if( maxcols == this%mincols() ) then
    if( size(this%access%values_) > 0 ) then
      call c_f_pointer (c_loc(this%access%values_(1)), data, [maxcols,int(this%access%rows_,c_int)])
      if( present(ncols) ) then
        ncols = maxcols
      endif
    endif
  else
    write(0,*) "ERROR: Cannot point connectivity pointer data(:,:) to irregular table"
  endif
end subroutine


subroutine atlas_Connectivity__row(this, row_idx, row, cols)
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), intent(in) :: row_idx
  integer(c_int), pointer, intent(out) :: row(:)
  integer(c_int), intent(out) ::  cols
  row  => this%access%values_(this%access%displs_(row_idx)+1:this%access%displs_(row_idx+1)+1)
  cols =  this%access%cols(row_idx)
end subroutine

subroutine atlas_Connectivity__add_values_args_size_t(this,rows,cols,values)
  use atlas_connectivity_c_binding
  class(atlas_Connectivity), intent(in) :: this
  integer(c_size_t) :: rows
  integer(c_size_t) :: cols
  integer(c_int) :: values(:)
  call atlas__connectivity__add_values(this%c_ptr(),rows,cols,values)
end subroutine

subroutine atlas_Connectivity__add_values_args_int(this,rows,cols,values)
  use atlas_connectivity_c_binding
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int) :: rows
  integer(c_int) :: cols
  integer(c_int) :: values(:)
  call atlas__connectivity__add_values(this%c_ptr(),int(rows,c_size_t),int(cols,c_size_t),values)
end subroutine

subroutine atlas_Connectivity__add_missing_args_size_t(this,rows,cols)
  use atlas_connectivity_c_binding
  class(atlas_Connectivity), intent(in) :: this
  integer(c_size_t) :: rows
  integer(c_size_t) :: cols
  call atlas__connectivity__add_missing(this%c_ptr(),rows,cols)
end subroutine

subroutine atlas_Connectivity__add_missing_args_int(this,rows,cols)
  use atlas_connectivity_c_binding
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int) :: rows
  integer(c_int) :: cols
  call atlas__connectivity__add_missing(this%c_ptr(),int(rows,c_size_t),int(cols,c_size_t))
end subroutine

!========================================================

function MultiBlockConnectivity_cptr(cptr) result(this)
  use atlas_connectivity_c_binding
  type(atlas_MultiBlockConnectivity) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
  call setup_access(this)
end function

function MultiBlockConnectivity_constructor() result(this)
  use atlas_connectivity_c_binding
  type(atlas_MultiBlockConnectivity) :: this
  this = MultiBlockConnectivity_cptr( atlas__MultiBlockConnectivity__create() )
  call this%return()
end function

function atlas_MultiBlockConnectivity__blocks(this) result(val)
  use atlas_connectivity_c_binding
  integer(c_size_t) :: val
  class(atlas_MultiBlockConnectivity), intent(in) :: this
  val = atlas__MultiBlockConnectivity__blocks(this%c_ptr())
end function

function atlas_MultiBlockConnectivity__block(this,block_idx) result(block)
  use atlas_connectivity_c_binding
  type(atlas_BlockConnectivity) :: block
  class(atlas_MultiBlockConnectivity), intent(in) :: this
  integer(c_size_t) :: block_idx
  call block%reset_c_ptr( atlas__MultiBlockConnectivity__block(this%c_ptr(),block_idx-1) )
end function

!========================================================

function BlockConnectivity_cptr(cptr) result(this)
  use atlas_connectivity_c_binding
  type(atlas_BlockConnectivity) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
end function

subroutine atlas_BlockConnectivity__delete(this)
  use atlas_connectivity_c_binding
  class(atlas_BlockConnectivity), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__BlockConnectivity__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
end subroutine

subroutine atlas_BlockConnectivity__copy(this,obj_in)
  class(atlas_BlockConnectivity), intent(inout) :: this
  class(atlas_RefCounted),   target, intent(in) :: obj_in
end subroutine

subroutine atlas_BlockConnectivity__data(this,data)
  use atlas_connectivity_c_binding
  class(atlas_BlockConnectivity), intent(in) :: this
  integer(c_int), pointer, intent(out) :: data(:,:)
  type(c_ptr) :: data_cptr
  integer(c_size_t) :: rows
  integer(c_size_t) :: cols
  call atlas__BlockConnectivity__data(this%c_ptr(),data_cptr,rows,cols)
  call c_f_pointer (data_cptr, data, [cols,rows])
end subroutine

function atlas_BlockConnectivity__rows(this) result(val)
  use atlas_connectivity_c_binding
  integer(c_size_t) :: val
  class(atlas_BlockConnectivity), intent(in) :: this
  val = atlas__BlockConnectivity__rows(this%c_ptr())
end function

function atlas_BlockConnectivity__cols(this) result(val)
  use atlas_connectivity_c_binding
  integer(c_size_t) :: val
  class(atlas_BlockConnectivity), intent(in) :: this
  val = atlas__BlockConnectivity__cols(this%c_ptr())
end function

function atlas_BlockConnectivity__missing_value(this) result(val)
  use atlas_connectivity_c_binding
  integer(c_int) :: val
  class(atlas_BlockConnectivity), intent(in) :: this
  val = atlas__BlockConnectivity__missing_value(this%c_ptr()) + 1
end function

!========================================================

subroutine setup_access(connectivity)
  class(atlas_Connectivity), intent(inout) :: connectivity
  type(c_ptr) :: ctxt_update
  type(c_ptr) :: ctxt_delete
  integer(c_int) :: have_ctxt_update
  integer(c_int) :: have_ctxt_delete
  have_ctxt_update = atlas__connectivity__ctxt_update(connectivity%c_ptr(),ctxt_update)
  if( have_ctxt_update == 1 ) then
    call c_f_pointer(ctxt_update,connectivity%access)
  else
    allocate( connectivity%access )
    call atlas__connectivity__register_update( connectivity%c_ptr(), c_funloc(update_access_c), c_loc(connectivity%access) )
    connectivity%access%connectivity_ptr = connectivity%c_ptr()
    call update_access(connectivity%access)
  endif
  have_ctxt_delete = atlas__connectivity__ctxt_delete(connectivity%c_ptr(),ctxt_delete)
  if( have_ctxt_update /= 1 ) then
    call atlas__connectivity__register_delete( connectivity%c_ptr(), c_funloc(delete_access_c), c_loc(connectivity%access) )
  endif
end subroutine

subroutine update_access_c(this_ptr) bind(c)
 type(c_ptr), value :: this_ptr
 type(atlas_ConnectivityAccess), pointer :: this
 call c_f_pointer(this_ptr,this)
 call update_access(this)
end subroutine

subroutine update_access(this)
  use atlas_connectivity_c_binding
  type(atlas_ConnectivityAccess), intent(inout) :: this
  integer :: jrow

  type(c_ptr) :: values_cptr
  type(c_ptr) :: displs_cptr
  type(c_ptr) :: counts_cptr
  integer(c_size_t) :: values_size
  integer(c_size_t) :: displs_size
  integer(c_size_t) :: counts_size

  this%missing_value_ = atlas__Connectivity__missing_value(this%connectivity_ptr)
  call atlas__Connectivity__values(this%connectivity_ptr,values_cptr,values_size)
  call atlas__Connectivity__displs(this%connectivity_ptr,displs_cptr,displs_size)
  call atlas__Connectivity__counts(this%connectivity_ptr,counts_cptr,counts_size)

  call c_f_pointer(values_cptr, this%values_, (/values_size/))
  call c_f_pointer(displs_cptr, this%displs_, (/displs_size/))
  call c_f_pointer(counts_cptr, this%cols,    (/counts_size/))
  this%rows_   = atlas__Connectivity__rows(this%connectivity_ptr)
  if( associated( this%row ) ) deallocate(this%row)
  allocate( this%row(this%rows_) )
  this%maxcols_ = 0
  this%mincols_ = huge(this%maxcols_)
  do jrow=1,this%rows_
    this%row(jrow)%col => this%values_(this%displs_(jrow)+1:this%displs_(jrow+1)+1)
    this%maxcols_ = max(this%maxcols_, this%cols(jrow) )
    this%mincols_ = min(this%mincols_, this%cols(jrow) )
  enddo
  if( associated( this%padded_) ) call update_padded(this)
end subroutine

subroutine update_padded(this)
  class(atlas_ConnectivityAccess), intent(inout) :: this
  integer(c_size_t) :: jrow, jcol
  if( associated(this%padded_) ) deallocate(this%padded_)
  allocate(this%padded_(this%maxcols_,this%rows()))
  this%padded_(:,:) = this%missing_value_
  do jrow=1,this%rows()
    do jcol=1,this%cols(jrow)
      this%padded_(jcol,jrow) = this%value(jcol,jrow)
    enddo
  enddo
end subroutine

subroutine delete_access_c(this_ptr) bind(c)
 type(c_ptr), value :: this_ptr
 type(atlas_ConnectivityAccess), pointer :: this
 call c_f_pointer(this_ptr,this)
 call delete_access(this)
end subroutine

subroutine delete_access(this)
  use atlas_connectivity_c_binding
  type(atlas_ConnectivityAccess), intent(inout) :: this
  integer :: jrow
  if( associated( this%row ) )    deallocate(this%row)
  if( associated( this%padded_) ) deallocate(this%padded_)
end subroutine


end module atlas_connectivity_module

