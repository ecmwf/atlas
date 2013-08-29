
module mpdata_module
  use state_module
  implicit none

  type, public :: Solver
    real :: dt_stability
    class(State), pointer :: state
    integer :: iter = 0
  contains
    procedure, pass :: init => Solver__init
    procedure, pass :: step => Solver__step
  end type Solver


  type, public, extends(Solver) :: MPDATA_Solver
    type(Field), pointer :: vol ! dual mesh volumes
    type(Field), pointer :: S ! dual mesh edge-normals
    type(Field), pointer :: grad_D ! gradient of depth
    type(Field), pointer :: R, Rex, Qadv, Q0, D0, Q, D, V
    type(Grid),  pointer :: g
    real, private :: dt

  contains
    procedure, pass :: init => MPDATA_Solver__init
    procedure, pass :: step => MPDATA_Solver__step
    procedure, private, pass :: compute_gradient => MPDATA_Solver__compute_gradient
    procedure, private, pass :: compute_Rn   => MPDATA_Solver__compute_Rn
    procedure, private, pass :: implicit_solve => MPDATA_Solver__implicit_solve
    procedure, private, pass :: compute_advective_velocities => MPDATA_Solver__compute_advective_velocities

  end type MPDATA_Solver


  type, public :: Model
    class(Solver), pointer :: solver
    class(State), pointer  :: state
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
      write(0,*) "iter =",self%solver%iter, &
         & "  time = ",self%state%time, &
         & "   dt = ",min( self%solver%dt_stability, tmax-t0 )
    end do
  end subroutine Model__solve_time_step


  subroutine Solver__init(self,g)
    class(Solver), intent(inout) :: self
    class(Grid), intent(in), target :: g
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


! ==========================================================================
! MPDATA Solver subroutines
! ==========================================================================

  subroutine MPDATA_Solver__init(self,g)
    class(MPDATA_Solver), intent(inout) :: self
    class(Grid), intent(in), target :: g
    class(FunctionSpace), pointer :: vertices
    class(FunctionSpace), pointer :: faces

    write(0,*) "MPDATA_Solver::init(g)"  
    self%g => g
    vertices => g%function_space("vertices")
    faces => g%function_space("faces")
    self%S => faces%field("dual_face_normal")
    self%vol => vertices%field("dual_volume")
    self%D => vertices%field("depth")
    self%Q => vertices%field("momentum")
    self%grad_D => vertices%add_vector_field("depth_gradient")
    self%R => vertices%add_vector_field("rhs_forcing")
    self%Rex => vertices%add_vector_field("rhs_forcing_exact_part")
    self%Qadv => vertices%add_vector_field("momentum_advected")
    self%D0 => vertices%add_vector_field("depth_prev_iter")
    self%Q0 => vertices%add_vector_field("momentum_prev_iter")
    self%V => vertices%add_vector_field("advective_velocity")

  end subroutine MPDATA_Solver__init

  subroutine MPDATA_Solver__step(self,tmax)
    class(MPDATA_Solver), intent(inout) :: self
    real, intent(in) :: tmax
    integer :: inode

    self%dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + self%dt
    
    ! Flow based on paper 'Szmelter, Smolarkiewicz; 
    !    An edge-based unstructured mesh discretisation in 
    !    geospherical framework, JCP, 2010'
    ! 
    ! This medium-high-level only involves loops over data points
    ! if (self%iter == 0) then
    !    0. First iteration... special case"
    !        0.1  Evaluate R^n = R_exp^n + R_impl^n"
    !        0.2  Evaluate advective velocities V^(n+1/2)"
    ! end if
    ! 1. Evaluate Q^n + 0.5 dt R^n"
    ! 2. Evaluate \hat(Q) = MPDATA( Q^n + 0.5 dt R^n , V^(n+1/2) , G )"
    ! 3. Evaluate R_exp = func( \hat(Q), grad( \hat(Q) ) )"
    ! 4. Evaluate \hat( \hat(Q) ) = \hat(Q) + 0.5 dt R_exp"
    ! 5. Compute Q^(n+1,m=1,3), with R_impl^(m-1) stored 
    !         from previous iteration"
    !     5.1 Evaluate Q^(n+1,m) = \hat( \hat(Q) ) + 0.5 dt R_impl^(m-1)"
    !     5.2 Evaluate R_impl^m = func( Q^(n+1,m) )"
    ! 6. Evaluate new advective velocities V^(n+1/2) = func(Q, Q0)
    !    V^(n+1/2) = 1/2 (3 V^n - V^(n-1) )
    ! 7. Store Q in Q0"

    ! Basically we need functions for:
    ! - R_exp(Q, grad(Q), geom)
    ! - R_impl(Q, geom)
    ! - grad(Q, geom)
    ! - V(Q,Q0)
    ! - MPDATA(Q,V,G)
    ! These low-level functions might need edge loops
    
    ! Initialise previous time step
    if (self%iter == 1) then
      self%Q0%array = self%Q%array
      self%D0%array = self%D%array
      call self%compute_Rn()
    end if
    
    call self%compute_advective_velocities()

    ! Backup solution
    self%Q0%array = self%Q%array
    self%D0%array = self%D%array

    ! Mimick MPDATA routines
    do inode=1,self%g%nb_nodes
      self%D%array(inode,1)    = self%D%array(inode,1)
      self%Qadv%array(inode,1) = self%Q%array(inode,1)
      self%Qadv%array(inode,2) = self%Q%array(inode,2)
    end do

    call self%implicit_solve()


    self%iter = self%iter + 1
  end subroutine MPDATA_Solver__step

  subroutine MPDATA_Solver__compute_advective_velocities(self)
    class(MPDATA_Solver), intent(in)    :: self
    real :: y, r, Qx, Qy, Q0x, Q0y, D, D0
    integer :: inode
    
    r = 6371.22e+03
   
    do inode=1,self%g%nb_nodes
      y     = self%g%nodes(inode,2)
      Qx    = self%Q%array(inode,1)
      Qy    = self%Q%array(inode,2)
      D     = self%D%array(inode,1)
      Q0x   = self%Q0%array(inode,1)
      Q0y   = self%Q0%array(inode,2)
      D0    = self%D0%array(inode,1)
      
      self%V%array(inode,1) = ( 1.5*Qx/D - 0.5*Q0x/D0 )*r
      self%V%array(inode,2) = ( 1.5*Qy/D - 0.5*Q0y/D0 )*r*cos(y)
      
    end do
  end subroutine MPDATA_Solver__compute_advective_velocities



  subroutine MPDATA_Solver__compute_Rn(self)
    class(MPDATA_Solver), intent(in)    :: self
    real :: f0, f, x, y, r, g, Qx, Qy, D, dDdx, sin_y, cos_y
    integer :: inode
    
    g = 9.80616
    r = 6371.22e+03
    f0 = 1.4584e-04 ! coriolis force (=2xearth's omega)

    call self%compute_gradient( self%D%array(:,1), self%grad_D%array )
    
    do inode=1,self%g%nb_nodes
      y     = self%g%nodes(inode,2)
      sin_y = sin(y)
      cos_y = cos(y)
      Qx    = self%Q%array(inode,1)
      Qy    = self%Q%array(inode,2)
      D     = self%D%array(inode,1)
      dDdx  = self%grad_D%array(inode,1)
      f     = f0 * sin_y
      
      self%Rex%array(inode,1) = -g/(r*cos_y)*D*dDdx
      self%Rex%array(inode,2) = -g/r*D*dDdx
      self%R%array(inode,1)   = self%Rex%array(inode,1) + f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
      self%R%array(inode,2)   = self%Rex%array(inode,2) - f*Qx - sin_y/(r*cos_y)*Qx*Qx/D
    end do
  end subroutine MPDATA_Solver__compute_Rn


  subroutine MPDATA_Solver__implicit_solve(self)
    class(MPDATA_Solver), intent(in)    :: self
    real :: f0, f, y, r, g, Qx, Qy, D, dDdx, sin_y, cos_y
    integer :: inode, m

    g = 9.80616
    r = 6371.22e+03
    f0 = 1.4584e-04 ! coriolis force (=2xearth's omega)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    call self%compute_gradient( self%D%array(:,1), self%grad_D%array )
    
    do inode=1,self%g%nb_nodes
      y  = self%g%nodes(inode,2)
      Qx = self%Q%array(inode,1)
      Qy = self%Q%array(inode,2)
      D  = self%D%array(inode,1)
      dDdx = self%grad_D%array(inode,1)
      self%Rex%array(inode,1) = -g/(r*cos(y))*D*dDdx
      self%Rex%array(inode,2) = -g/r*D*dDdx
    end do

    do inode=1,self%g%nb_nodes
      y     = self%g%nodes(inode,2)
      sin_y = sin(y)
      cos_y = cos(y)
      D = self%D%array(inode,1)
      f = f0 * sin_y

      do m=1,3 ! Three iterations at most is enough to converge
        Qx = self%Q%array(inode,1)
        Qy = self%Q%array(inode,2)

        self%R%array(inode,1) = self%Rex%array(inode,1) + f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
        self%R%array(inode,2) = self%Rex%array(inode,2) - f*Qx - sin_y/(r*cos_y)*Qx*Qx/D
      
        self%Q%array(inode,1) = self%Qadv%array(inode,1) + 0.5*self%dt*self%R%array(inode,1)
        self%Q%array(inode,2) = self%Qadv%array(inode,2) + 0.5*self%dt*self%R%array(inode,2)
      end do
    
    end do
  end subroutine MPDATA_Solver__implicit_solve

!================================================================
! Subroutines independent of equations used
!================================================================
  subroutine MPDATA_Solver__compute_gradient(self,Q,gradQ)
    class(MPDATA_Solver), intent(in)    :: self
    real, dimension(:),   intent(in)    :: Q
    real, dimension(:,:), intent(inout) :: gradQ
    
    real    :: sx,sy,avgQ
    integer :: ip,ied,p1,p2, nb_internal_edges, nb_edges

    nb_internal_edges = self%g%nb_faces
    nb_edges = self%g%nb_faces

    ! derivatives 
    do ip = 1,size(Q)
      gradQ(ip,1) = 0.
      gradQ(ip,2) = 0.
    enddo

    do ied = 1,nb_internal_edges
      p1  = self%g%faces(ied,1)
      p2  = self%g%faces(ied,2)
      sx  = self%S%array(ied,1)
      sy  = self%S%array(ied,2)
      avgQ = 0.5*( Q(p1) + Q(p2) )
      gradQ(p1,1) = gradQ(p1,1) + sx*avgQ
      gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
      gradQ(p2,1) = gradQ(p2,1) - sx*avgQ
      gradQ(p2,2) = gradQ(p2,2) - sy*avgQ
    enddo

    ! special treatment for the north & south pole cell faces 
    do ied = nb_internal_edges+1,nb_edges
      p1  = self%g%faces(ied,1)
      p2  = self%g%faces(ied,2)
      sx  = self%S%array(ied,1)
      sy  = self%S%array(ied,2)
      avgQ = 0.5*( Q(p1) + Q(p2) )
      gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
      gradQ(p2,2) = gradQ(p2,2) + sy*avgQ
    enddo

  end subroutine MPDATA_Solver__compute_gradient

end module mpdata_module

