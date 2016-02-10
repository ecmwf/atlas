
module atlas_connectivity_module

use iso_c_binding, only : c_funptr, c_ptr, c_loc, c_f_pointer, c_f_procpointer, c_funloc, c_int, c_size_t
use atlas_refcounted_module
implicit none

private :: c_funptr
private :: c_ptr
private :: c_loc
private :: c_f_pointer
private :: c_funloc
private :: c_int
private :: c_size_t

private

!-----------------------------
! atlas_Connectivity         !
!-----------------------------

type, extends(atlas_refcounted), public :: atlas_Connectivity

! Public members
  type( atlas_ConnectivityAccess ), public, pointer :: access

contains

! Public methods
  procedure, public :: value    => atlas_Connectivity__value
  procedure, public :: rows     => atlas_Connectivity__rows
  procedure, public :: cols     => atlas_Connectivity__cols
  procedure, public :: maxcols  => atlas_Connectivity__maxcols
  procedure, public :: mincols  => atlas_Connectivity__mincols
  procedure, public :: data     => atlas_Connectivity__data

  procedure, public :: delete => atlas_Connectivity__delete

end type

interface atlas_Connectivity
  module procedure Connectivity_ctr
end interface

!-------------------------------
! Helper types                 !
!-------------------------------

type, public :: atlas_ConnectivityAccessRow
  integer, pointer :: col(:)
contains
end type

type, public :: atlas_ConnectivityAccess
  integer(c_int), private, pointer :: values_(:) => null()
  integer(c_size_t), private, pointer :: displs_(:) => null()
  integer(c_size_t), public,  pointer :: cols(:)    => null()
  type(atlas_ConnectivityAccessRow), public, pointer :: row(:) => null()
  integer(c_size_t), private :: rows_
  integer, private, pointer :: padded_(:,:) => null()
  integer(c_size_t), private :: maxcols_, mincols_
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

function atlas_Connectivity__cptr(cptr) result(connectivity)
  type(atlas_Connectivity) :: connectivity
  type(c_ptr), intent(in) :: cptr
  call connectivity%reset_c_ptr( cptr )
end function

function Connectivity_ctr() result(connectivity)
  use atlas_connectivity_c_binding
  type(atlas_Connectivity) :: connectivity
  type(c_ptr) :: ctxt_update
  type(c_ptr) :: ctxt_delete
  integer(c_int) :: have_ctxt_update
  integer(c_int) :: have_ctxt_delete
  connectivity = atlas_Connectivity__cptr( atlas__Connectivity__create() )

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


  call connectivity%return()
end function

subroutine atlas_Connectivity__delete(this)
  use atlas_connectivity_c_binding
  class(atlas_Connectivity), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Connectivity__delete(this%c_ptr())
  end if
  call this%reset_c_ptr()
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

pure function atlas_Connectivity__value(this,c,r) result(val)
  integer(c_int) :: val
  class(atlas_Connectivity), intent(in) :: this
  integer(c_size_t), intent(in) :: r,c
  val = this%access%values_(c+this%access%displs_(r))
end function

pure function atlas_Connectivity__rows(this) result(val)
  integer(c_size_t) :: val
  class(atlas_Connectivity), intent(in) :: this
  val = this%access%rows_
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

subroutine atlas_Connectivity__data(this, padded, cols)
  class(atlas_Connectivity), intent(in) :: this
  integer(c_int), pointer, intent(out) :: padded(:,:)
  integer(c_size_t), pointer, intent(out), optional :: cols(:)
  if( .not. associated(this%access%padded_) ) call update_padded(this%access)
  padded => this%access%padded_
  if( present(cols) ) cols => this%access%cols
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

  call atlas__Connectivity__values(this%connectivity_ptr,values_cptr,values_size)
  call atlas__Connectivity__displs(this%connectivity_ptr,displs_cptr,displs_size)
  call atlas__Connectivity__counts(this%connectivity_ptr,counts_cptr,counts_size)

  call c_f_pointer(values_cptr, this%values_, (/values_size/))
  call c_f_pointer(displs_cptr, this%displs_, (/displs_size/))
  call c_f_pointer(counts_cptr, this%cols,    (/counts_size/))
  this%rows_   = counts_size

  if( associated( this%row ) ) deallocate(this%row)
  allocate( this%row(this%rows()) )
  this%maxcols_ = 0
  this%mincols_ = huge(this%maxcols_)
  do jrow=1,this%rows()
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

