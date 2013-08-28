
module mpdata_module
  use state_module
  implicit none

  type, public :: Solver
    real :: dt_stability
    class(State), pointer :: state
  contains
    procedure, pass :: step => Solver__step
  end type Solver


  type, public, extends(Solver) :: MPDATA_Solver
  contains
  end type MPDATA_Solver


  type, public :: Model
    class(Solver), pointer :: solver
    class(State), pointer  :: state
    integer :: iter = 0
  contains
    procedure, pass :: solve_time_step => Model__solve_time_step
  end type Model

contains

  subroutine Model__solve_time_step(self,dt)
    class(Model), intent(inout) :: self
    real, intent(in) :: dt
    real :: tmax, t0

    tmax = self%state%time+dt
    do while (self%state%time < tmax)
      t0 = self%state%time
      call self%solver%step(tmax)
      self%iter = self%iter + 1
      write(0,*) "iter =",self%iter,"  time = ",self%state%time,"   dt = ",min( self%solver%dt_stability, tmax-t0 )
    end do
  end subroutine Model__solve_time_step

  subroutine Solver__step(self,tmax)
    class(Solver), intent(inout) :: self
    real, intent(in) :: tmax
    real :: dt
    dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + dt
  end subroutine Solver__step
end module mpdata_module

