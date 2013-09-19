! State module
! -----------------
! This module contains the State class which contains
! functionspaces and fields, time-level, ...
module model_module
  use grid_module
  implicit none
  public
 
  type, public :: State_class
    character(len=30) :: name
    real(kind=jprb)              :: time
    class(FunctionSpace_class), pointer :: function_space
    type(FieldPtr), dimension(:), allocatable :: fields

  contains
    procedure, pass :: init => State__init
    procedure, pass :: destruct => State__destruct
    procedure, pass :: add_field => State__add_field
    procedure, pass :: field => State__field
    procedure, pass :: has_field => State__has_field
  end type State_class
  
  type, public :: StatePtr
    type(State_class), pointer :: ptr
  end type StatePtr
  

 type, public :: Solver_class
    real(kind=jprb) :: dt_stability
    class(State_class), pointer :: state
    class(Model_class), pointer :: model
    integer :: iter = 0
  contains
    procedure, pass :: init => Solver__init
    procedure, pass :: step => Solver__step
  end type Solver_class


  type, public :: Model_class
    class(Solver_class), pointer :: solver
    class(State_class),  pointer :: state
    class(Grid_class),   pointer :: grid
    class(Model_class),  pointer :: ptr 

  contains
    procedure, pass :: init => Model__init
    procedure, pass :: solve_time_step => Model__solve_time_step
  end type Model_class

contains

! ------------------------------------------------------------------------------------
!                                  State subroutines
!-------------------------------------------------------------------------------------

  subroutine State__init(self, name, function_space)
    class(State_class), intent(inout) :: self
    class(FunctionSpace_class), intent(in), target ::function_space
    character(len=*), intent(in) :: name
    write(0,*) "State::init(",name,")"  
    self%name = name  
    self%function_space => function_space
    allocate(self%fields(0))
  end subroutine State__init
  
  subroutine State__destruct(self)
    class(State_class), intent(inout) :: self
    write(0,*) "State::destruct"
    deallocate(self%fields)
  end subroutine State__destruct
  
  subroutine State__add_field(self,field)
    class(State_class), intent(inout)   :: self
    class(Field_class), pointer :: field
    type(FieldPtr), allocatable :: tmp(:)
    call move_alloc(self%fields,tmp)
    allocate(self%fields(size(tmp)+1))
    self%fields(:size(tmp)) = tmp
    self%fields(size(self%fields))%ptr => field
    deallocate(tmp)
  end subroutine State__add_field

  function State__field(self, name) result(field)
    class(State_class), intent(in) :: self
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: field
    integer :: f
    do f=1,size(self%fields)
      field => self%fields(f)%ptr
      if( field%name == name) then
        return
      end if
    end do
    write(0,*) 'Could not find field ',name, ' in ',self%name
    call abort
  end function State__field

  logical function State__has_field(self, name)
    class(State_class), intent(in) :: self
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: field
    logical :: has_field
    integer :: f
    do f=1,size(self%fields)
      field => self%fields(f)%ptr
      if( field%name == name) then
        State__has_field = .True.
        return
      end if
    end do
    State__has_field = .False.
  end function State__has_field

! ------------------------------------------------------------------------------------
!                                  Model subroutines
!-------------------------------------------------------------------------------------

  subroutine Model__init(self, grid)
    class(Model_class), intent(inout), target :: self
    class(Grid_class), intent(in), target :: grid
    self%ptr => self
    self%grid => grid
  end subroutine Model__init

  subroutine Model__solve_time_step(self,dt)
    class(Model_class), intent(inout) :: self
    real(kind=jprb), intent(in) :: dt
    real(kind=jprb) :: tmax, t0

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

  subroutine Solver__init(self,model)
    class(Solver_class), intent(inout) :: self
    class(Model_class), intent(in), target :: model
    self%model => model
    write(0,*) "Solver::init(g)"  
  end subroutine Solver__init


  subroutine Solver__step(self,tmax)
    class(Solver_class), intent(inout) :: self
    real(kind=jprb), intent(in) :: tmax
    real(kind=jprb) :: dt
    dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + dt
    self%iter = self%iter + 1
  end subroutine Solver__step


end module model_module