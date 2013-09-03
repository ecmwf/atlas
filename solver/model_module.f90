! State module
! -----------------
! This module contains the State class which contains
! functionspaces and fields, time-level, ...
module model_module
  use grid_module
  implicit none
  public
 
  type, public :: State
    character(len=30) :: name
    real              :: time
    class(FunctionSpace), pointer :: function_space
    type(FieldPtr), dimension(:), allocatable :: fields

  contains
    procedure, pass :: init => State__init
    procedure, pass :: destruct => State__destruct
    procedure, pass :: add_field => State__add_field
    procedure, pass :: field => State__field
  end type State
  
  type, public :: StatePtr
    type(State), pointer :: ptr
  end type StatePtr
  

 type, public :: Solver
    real :: dt_stability
    class(State), pointer :: state
    class(Model), pointer :: model
    integer :: iter = 0
  contains
    procedure, pass :: init => Solver__init
    procedure, pass :: step => Solver__step
  end type Solver


  type, public :: Model
    class(Solver), pointer :: solver
    class(State),  pointer :: state
    class(Grid),   pointer :: grid
    class(Model),  pointer :: ptr 

  contains
    procedure, pass :: init => Model__init
    procedure, pass :: solve_time_step => Model__solve_time_step
  end type Model

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

  function State__field(self, name) result(field_)
    class(State), intent(in) :: self
    character(len=*), intent(in) :: name
    type(Field), pointer :: field_
    integer :: f
    do f=1,size(self%fields)
      field_ => self%fields(f)%ptr
      if( field_%name == name) then
        return
      end if
    end do
    !abort("No field named "//trim(name)//" in function_space")
  end function State__field

! ------------------------------------------------------------------------------------
!                                  Model subroutines
!-------------------------------------------------------------------------------------

  subroutine Model__init(self, g)
    class(Model), intent(inout), target :: self
    class(Grid), intent(in), target :: g
    self%ptr => self
    self%grid => g
  end subroutine Model__init

  subroutine Model__solve_time_step(self,dt)
    class(Model), intent(inout) :: self
    real, intent(in) :: dt
    real :: tmax, t0

    tmax = self%state%time+dt
    do while (self%state%time < tmax)
      t0 = self%state%time
      call self%solver%step(tmax)
      write(0,*) "iter =",self%solver%iter, &
         & "  time = ",self%state%time, &
         & "   dt = ",min( self%solver%dt_stability, tmax-t0 )
    end do
  end subroutine Model__solve_time_step

! ------------------------------------------------------------------------------------
!                                  Solver subroutines
!-------------------------------------------------------------------------------------

  subroutine Solver__init(self,model_)
    class(Solver), intent(inout) :: self
    class(Model), intent(in), target :: model_
    self%model => model_
    write(0,*) "Solver::init(g)"  
  end subroutine Solver__init


  subroutine Solver__step(self,tmax)
    class(Solver), intent(inout) :: self
    real, intent(in) :: tmax
    real :: dt
    dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + dt
    self%iter = self%iter + 1
  end subroutine Solver__step


end module model_module