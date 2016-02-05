
module atlas_connectivity_module

implicit none

type, public :: atlas_ConnectivityRow
  integer, pointer :: col(:)
contains
end type

type, public :: atlas_Connectivity
  integer, pointer :: values_(:)
  integer, pointer :: displs_(:)
  integer, pointer :: cols(:)
  type(atlas_ConnectivityRow), allocatable :: row(:)
  integer :: rows
contains
  procedure, public :: value
  procedure, private :: update
end type

interface atlas_Connectivity
  module procedure Connectivity_ctr
end interface

contains

function Connectivity_ctr(values,displs,counts) result(conn)
  type(atlas_Connectivity) :: conn
  integer, target :: values(:)
  integer, target :: displs(:)
  integer, target :: counts(:)
  conn%values_ => values
  conn%displs_ => displs
  conn%cols    => counts
  conn%rows = size(displs-1)
  call conn%update()
end function

pure function value(this,c,r) result(val)
  integer :: val
  class(atlas_Connectivity), intent(in) :: this
  integer, intent(in) :: r,c 
  val = this%values_(c+this%displs_(r))
end function

subroutine update(this)
  class(atlas_Connectivity) :: this
  integer :: jrow
  allocate(this%row(this%rows))
  if( allocated( this%row ) ) deallocate(this%row)
  allocate( this%row(this%rows) )
  do jrow=1,this%rows
    this%row(jrow)%col => this%values_(this%displs_(jrow)+1:this%displs_(jrow+1)+1)
  enddo
end subroutine

end module atlas_connectivity_module
