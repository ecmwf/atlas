! State module
! -----------------
! This module contains the State class which contains
! functionspaces and fields, time-level, ...
module state_module
  use grid_module
  implicit none
  public
  
  type, public :: State
    character(len=30) :: name
    real              :: time
    type(FunctionSpace), pointer :: function_space
    type(FieldPtr), dimension(:), allocatable :: fields

  contains
    procedure, pass :: init => State__init
    procedure, pass :: destruct => State__destruct
    procedure, pass :: add_field => State__add_field
  end type State
  
  type, public :: StatePtr
    type(State), pointer :: ptr
  end type StatePtr
  
  contains
  
! ------------------------------------------------------------------------------------
!                                  State subroutines
!-------------------------------------------------------------------------------------

  subroutine State__init(self, name, function_space)
    class(State), intent(inout) :: self
    class(FunctionSpace), intent(in), target ::function_space
    character(len=*), intent(in) :: name
    write(0,*) "State::init(",name,")"  
    self%name = name  
    self%function_space => function_space
    allocate(self%fields(0))
  end subroutine State__init
  
  subroutine State__destruct(self)
    class(State), intent(inout) :: self
    write(0,*) "State::destruct"
    deallocate(self%fields)
  end subroutine State__destruct
  
  subroutine State__add_field(self,field_)
    class(State), intent(inout)   :: self
    type(Field), pointer :: field_
    type(FieldPtr), allocatable :: tmp(:)
    call move_alloc(self%fields,tmp)
    allocate(self%fields(size(tmp)+1))
    self%fields(:size(tmp)) = tmp
    self%fields(size(self%fields))%ptr => field_
    deallocate(tmp)
  end subroutine State__add_field

end module state_module