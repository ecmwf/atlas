
module mpdata_module
  use model_module
  implicit none

  type, public, extends(Solver) :: MPDATA_Solver
    type(Field), pointer :: vol ! dual mesh volumes
    type(Field), pointer :: S ! dual mesh edge-normals
    real :: dt

  contains
    procedure, public,  pass :: init => MPDATA_Solver__init
    procedure, public,  pass :: step => MPDATA_Solver__step
    procedure, pass :: compute_gradient => MPDATA_Solver__compute_gradient
    procedure, pass :: compute_Rn   => MPDATA_Solver__compute_Rn
    procedure, pass :: implicit_solve => MPDATA_Solver__implicit_solve
    procedure, pass :: compute_advective_velocities => MPDATA_Solver__compute_advective_velocities
    procedure, pass :: mpdata_gage => MPDATA_Solver__mpdata_gage
    procedure, pass :: backup_solution => MPDATA_Solver__backup_solution
    procedure, pass :: add_forcing_to_solution => MPDATA_Solver__add_forcing_to_solution
    procedure, pass :: advect_solution => MPDATA_Solver__advect_solution
    
  end type MPDATA_Solver


contains


! ==========================================================================
! MPDATA Solver subroutines
! ==========================================================================

  function new_MPDATA_Solver(model_) result(solver_)
    class(Model), pointer :: model_
    class(Solver), pointer :: solver_
    allocate( MPDATA_Solver :: solver_ )
    call solver_%init(model_)
  end function new_MPDATA_Solver

  subroutine MPDATA_Solver__init(self,model_)
    class(MPDATA_Solver), intent(inout) :: self
    class(Model), intent(in), target :: model_
    class(FunctionSpace), pointer :: vertices
    class(FunctionSpace), pointer :: faces

    call self%Solver%init(model_)
    write(0,*) "MPDATA_Solver::init(g)"  
    vertices => self%model%grid%function_space("vertices")
    faces => self%model%grid%function_space("faces")
    self%S => faces%field("dual_face_normal")
    self%vol => vertices%field("dual_volume")
  end subroutine MPDATA_Solver__init


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
  subroutine MPDATA_Solver__step(self,tmax)
    class(MPDATA_Solver), intent(inout) :: self
    real, intent(in) :: tmax
    integer :: inode

    self%dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + self%dt
    
    call self%backup_solution()
    
    ! Initialise previous time step
    if (self%iter == 0) then
      call self%compute_Rn()
    end if
    
    call self%compute_advective_velocities()
    call self%add_forcing_to_solution()
    call self%advect_solution()
    call self%implicit_solve()
    self%iter = self%iter + 1
  end subroutine MPDATA_Solver__step

  !================================================================
  ! Subroutines that need to be implemented
  !================================================================

  subroutine MPDATA_Solver__backup_solution(self)
    class(MPDATA_Solver), intent(in) :: self
  end subroutine MPDATA_Solver__backup_solution

  subroutine MPDATA_Solver__add_forcing_to_solution(self)
    class(MPDATA_Solver), intent(in) :: self
  end subroutine MPDATA_Solver__add_forcing_to_solution

  subroutine MPDATA_Solver__compute_advective_velocities(self)
    class(MPDATA_Solver), intent(in)    :: self
  end subroutine MPDATA_Solver__compute_advective_velocities

  subroutine MPDATA_Solver__advect_solution(self)
    class(MPDATA_Solver), intent(in)    :: self
  end subroutine MPDATA_Solver__advect_solution

  subroutine MPDATA_Solver__compute_Rn(self)
    class(MPDATA_Solver), intent(in)    :: self
  end subroutine MPDATA_Solver__compute_Rn

  subroutine MPDATA_Solver__implicit_solve(self)
    class(MPDATA_Solver), intent(in)    :: self
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

    nb_internal_edges = self%model%grid%nb_faces
    nb_edges = self%model%grid%nb_faces

    ! derivatives 
    do ip = 1,size(Q)
      gradQ(ip,1) = 0.
      gradQ(ip,2) = 0.
    enddo

    do ied = 1,nb_internal_edges
      p1  = self%model%grid%faces(ied,1)
      p2  = self%model%grid%faces(ied,2)
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
      p1  = self%model%grid%faces(ied,1)
      p2  = self%model%grid%faces(ied,2)
      sx  = self%S%array(ied,1)
      sy  = self%S%array(ied,2)
      avgQ = 0.5*( Q(p1) + Q(p2) )
      gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
      gradQ(p2,2) = gradQ(p2,2) + sy*avgQ
    enddo

  end subroutine MPDATA_Solver__compute_gradient

  subroutine MPDATA_Solver__mpdata_gage(self,Q,V, Qadv)
    class(MPDATA_Solver), intent(in)  :: self
    real, dimension(:),   intent(inout)  :: Q
    real, dimension(:,:), intent(in)  :: V
    real, dimension(:),   intent(out) :: Qadv
    
    integer :: inode, nb_internal_edges, ied, pass, p1,p2, apos, aneg, nb_nodes
    real :: rhs(self%model%grid%nb_nodes)
    real :: aun(self%model%grid%nb_faces)
    real :: sx, sy, flux

    write(0,*) "mpdata"

    nb_internal_edges = self%model%grid%nb_faces
    nb_nodes = self%model%grid%nb_nodes

    do pass=1,1 ! 2 or more passes here!!!
      
      do ied=1,nb_internal_edges
        p1 = self%model%grid%faces(ied,1)
        p2 = self%model%grid%faces(ied,2)
        sx = self%S%array(ied,1)
        sy = self%S%array(ied,2)

        aun(ied) = 0.5*(V(p1,1)+V(p2,1))*sx + 0.5*(V(p1,2)+V(p2,2))*sy
        apos = max(0.,aun(ied))
        aneg = min(0.,aun(ied))
        flux = Q(p1)*apos + Q(p2)*aneg
        rhs(p1) = rhs(p1) + flux
        rhs(p2) = rhs(p2) - flux
      end do

      do inode=1,nb_nodes
        Q(inode) = Q(inode) - rhs(inode) * self%dt / self%vol%array(inode,1)
      end do
    end do
    
    Qadv = Q
    
  end subroutine MPDATA_Solver__mpdata_gage

end module mpdata_module

