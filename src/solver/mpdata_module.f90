
module mpdata_module
  use model_module
  implicit none

  type, public, extends(Solver_class) :: MPDATA_Solver
    class(Field_class), pointer :: vol ! dual mesh volumes
    class(Field_class), pointer :: S ! dual mesh edge-normals
    real, dimension(:,:), allocatable :: gradQ
    real :: dt

  contains
    procedure, public,  pass :: init => MPDATA_Solver__init
    procedure, public,  pass :: step => MPDATA_Solver__step
    procedure, pass :: compute_gradient => MPDATA_Solver__compute_gradient
    procedure, pass :: compute_Rn   => MPDATA_Solver__compute_Rn
    procedure, pass :: implicit_solve => MPDATA_Solver__implicit_solve
    procedure, pass :: compute_advective_velocities => MPDATA_Solver__compute_advective_velocities
    procedure, pass :: mpdata_gauge => MPDATA_Solver__mpdata_gauge
    procedure, pass :: mpdata_abs => MPDATA_Solver__mpdata_abs
    procedure, pass :: backup_solution => MPDATA_Solver__backup_solution
    procedure, pass :: add_forcing_to_solution => MPDATA_Solver__add_forcing_to_solution
    procedure, pass :: advect_solution => MPDATA_Solver__advect_solution
    procedure, pass :: compute_weighted_average => MPDATA_Solver__compute_weighted_average

    
  end type MPDATA_Solver


contains


! ==========================================================================
! MPDATA Solver subroutines
! ==========================================================================

  function new_MPDATA_Solver(model) result(solver)
    class(Model_class), pointer :: model
    class(Solver_class), pointer :: solver
    allocate( MPDATA_Solver :: solver )
    call solver%init(model)
  end function new_MPDATA_Solver

  subroutine MPDATA_Solver__init(self,model)
    class(MPDATA_Solver), intent(inout) :: self
    class(Model_class), intent(in), target :: model
    class(FunctionSpace_class), pointer :: vertices
    class(FunctionSpace_class), pointer :: faces

    call self%Solver_class%init(model)
    write(0,*) "MPDATA_Solver::init(g)"  
    vertices => self%model%grid%function_space("vertices")
    faces => self%model%grid%function_space("faces")
    self%S => faces%field("dual_face_normal")
    self%vol => vertices%field("dual_volume")
    allocate( self%gradQ(vertices%nb_nodes,self%model%grid%dimension) )
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
    class(MPDATA_Solver), intent(inout) :: self
  end subroutine MPDATA_Solver__backup_solution

  subroutine MPDATA_Solver__add_forcing_to_solution(self)
    class(MPDATA_Solver), intent(inout) :: self
  end subroutine MPDATA_Solver__add_forcing_to_solution

  subroutine MPDATA_Solver__compute_advective_velocities(self)
    class(MPDATA_Solver), intent(inout)    :: self
  end subroutine MPDATA_Solver__compute_advective_velocities

  subroutine MPDATA_Solver__advect_solution(self)
    class(MPDATA_Solver), intent(inout)    :: self
  end subroutine MPDATA_Solver__advect_solution

  subroutine MPDATA_Solver__compute_Rn(self)
    class(MPDATA_Solver), intent(inout)    :: self
  end subroutine MPDATA_Solver__compute_Rn

  subroutine MPDATA_Solver__implicit_solve(self)
    class(MPDATA_Solver), intent(inout)    :: self
  end subroutine MPDATA_Solver__implicit_solve

!================================================================
! Subroutines independent of equations used
!================================================================

  subroutine MPDATA_Solver__compute_gradient(self,Q,gradQ)
    class(MPDATA_Solver), intent(inout) :: self
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
    end do

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
    end do

    ! special treatment for the north & south pole cell faces 
    do ied = nb_internal_edges+1,nb_edges
      p1  = self%model%grid%faces(ied,1)
      p2  = self%model%grid%faces(ied,2)
      sx  = self%S%array(ied,1)
      sy  = self%S%array(ied,2)
      avgQ = 0.5*( Q(p1) + Q(p2) )
      gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
      gradQ(p2,2) = gradQ(p2,2) + sy*avgQ
    end do

  end subroutine MPDATA_Solver__compute_gradient

  subroutine MPDATA_Solver__compute_weighted_average(self,Q,avgQ,absS)
    class(MPDATA_Solver), intent(inout) :: self
    real, dimension(:),   intent(in)    :: Q
    real, dimension(:),   intent(inout) :: avgQ
    real, dimension(:),   intent(inout) :: absS

    real    :: sx,sy,Qf, normS
    integer :: ip,ied,p1,p2, nb_internal_edges, nb_edges

    nb_internal_edges = self%model%grid%nb_faces
    nb_edges = self%model%grid%nb_faces

    ! derivatives 
    avgQ = 0
    absS = 0

    do ied = 1,nb_internal_edges
      p1  = self%model%grid%faces(ied,1)
      p2  = self%model%grid%faces(ied,2)
      sx  = self%S%array(ied,1)
      sy  = self%S%array(ied,2)
      Qf = 0.5*( Q(p1) + Q(p2) )
      normS = sqrt(sx*sx+sy*sy)
      avgQ(p1) = avgQ(p1) + normS*Qf
      avgQ(p2) = avgQ(p2) + normS*Qf
      absS(p1) = absS(p1) + normS
      absS(p2) = absS(p2) + normS
    end do
  end subroutine MPDATA_Solver__compute_weighted_average


  subroutine MPDATA_Solver__mpdata_gauge(self,Q,V)
    class(MPDATA_Solver), intent(inout)  :: self
    real, dimension(:),   intent(inout)  :: Q
    real, dimension(:,:), intent(in)  :: V
    
    integer :: inode, nb_internal_edges, ied, pass, p1,p2, apos, aneg, nb_nodes
    real :: advection(self%model%grid%nb_nodes)
    real :: aun(self%model%grid%nb_faces)
    real :: Qabs(self%model%grid%nb_nodes)
    real :: Sabs(self%model%grid%nb_nodes)
    real :: Qavg(self%model%grid%nb_nodes)

    real :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, eps
    real :: denom
    class(Field_class), pointer :: Dgrad

    eps = 1e-10
    nb_internal_edges = self%model%grid%nb_faces
    nb_nodes = self%model%grid%nb_nodes

    ! 1. First pass
    ! -------------
      
    ! Compute the rhs in vertices, and normal velocity in faces
    advection = 0.
    do ied=1,nb_internal_edges
      p1 = self%model%grid%faces(ied,1)
      p2 = self%model%grid%faces(ied,2)
      sx = self%S%array(ied,1)
      sy = self%S%array(ied,2)

      aun(ied) = 0.5*(V(p1,1)+V(p2,1))*sx + 0.5*(V(p1,2)+V(p2,2))*sy
      apos = max(0.,aun(ied))
      aneg = min(0.,aun(ied))
      flux = Q(p1)*apos + Q(p2)*aneg
      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode) * self%dt / self%vol%array(inode,1)
      Qabs(inode) = abs(Q(inode))
    end do

  
    ! 2. Second pass
    ! -------------
    
    ! Compute derivatives for mpdata
    call self%compute_gradient(Qabs,self%gradQ)

    ! Compute antidiffusive normal velocity in faces
    advection = 0.
    do ied=1,nb_internal_edges
      p1 = self%model%grid%faces(ied,1)
      p2 = self%model%grid%faces(ied,2)
      sx = self%S%array(ied,1)
      sy = self%S%array(ied,2)

      ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
      volume_of_two_cells = self%vol%array(p1,1) + self%vol%array(p2,1)
      dQdx = (self%gradQ(p1,1)+self%gradQ(p2,1)) / volume_of_two_cells
      dQdy = (self%gradQ(p1,2)+self%gradQ(p2,2)) / volume_of_two_cells
      Vx = 0.5*(V(p1,1)+V(p2,1))
      Vy = 0.5*(V(p1,2)+V(p2,2))

      ! variable sign option with asymptotic analysis, (mpdata gauge)
      aun(ied) = abs(aun(ied))*(Q(p2)-Q(p1))*0.5 - 0.5*self%dt*aun(ied)*(Vx*dQdx+Vy*dQdy)
      flux = aun(ied)

      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode) * self%dt / self%vol%array(inode,1)
    end do

  end subroutine MPDATA_Solver__mpdata_gauge


  subroutine MPDATA_Solver__mpdata_abs(self,Q,V)
    class(MPDATA_Solver), intent(inout)  :: self
    real, dimension(:),   intent(inout)  :: Q
    real, dimension(:,:), intent(in)  :: V
    
    integer :: inode, nb_internal_edges, ied, pass, p1,p2, apos, aneg, nb_nodes
    real :: advection(self%model%grid%nb_nodes)
    real :: aun(self%model%grid%nb_faces)
    real :: Qabs(self%model%grid%nb_nodes)
    real :: Sabs(self%model%grid%nb_nodes)
    real :: Qavg(self%model%grid%nb_nodes)

    real :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, eps
    real :: denom
    class(Field_class), pointer :: Dgrad

    eps = 1e-10
    nb_internal_edges = self%model%grid%nb_faces
    nb_nodes = self%model%grid%nb_nodes

    ! 1. First pass
    ! -------------
      
    ! Compute the rhs in vertices, and normal velocity in faces
    advection = 0.
    do ied=1,nb_internal_edges
      p1 = self%model%grid%faces(ied,1)
      p2 = self%model%grid%faces(ied,2)
      sx = self%S%array(ied,1)
      sy = self%S%array(ied,2)

      aun(ied) = 0.5*(V(p1,1)+V(p2,1))*sx + 0.5*(V(p1,2)+V(p2,2))*sy
      apos = max(0.,aun(ied))
      aneg = min(0.,aun(ied))
      flux = Q(p1)*apos + Q(p2)*aneg
      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode) * self%dt / self%vol%array(inode,1)
      Qabs(inode) = abs(Q(inode))
    end do

  
    ! 2. Second pass
    ! -------------
    
    ! Compute derivatives for mpdata
    call self%compute_gradient(Qabs,self%gradQ)
    call self%compute_weighted_average(Q,Qavg,Sabs)


    ! Compute antidiffusive normal velocity in faces
    advection = 0.
    do ied=1,nb_internal_edges
      p1 = self%model%grid%faces(ied,1)
      p2 = self%model%grid%faces(ied,2)
      sx = self%S%array(ied,1)
      sy = self%S%array(ied,2)
      volume_of_two_cells = self%vol%array(p1,1) + self%vol%array(p2,1)
      dQdx = (self%gradQ(p1,1)+self%gradQ(p2,1)) / volume_of_two_cells
      dQdy = (self%gradQ(p1,2)+self%gradQ(p2,2)) / volume_of_two_cells
      Vx = 0.5*(V(p1,1)+V(p2,1))
      Vy = 0.5*(V(p1,2)+V(p2,2))

      ! variable sign option with absolute values
      denom = ( Qavg(p1)+Qavg(p2) ) / ( Sabs(p1) + Sabs(p2) )
      aun(ied) = abs(aun(ied))*( Qabs(p2)-Qabs(p1) )/( Qabs(p2)+Qabs(p1)+eps ) - 0.5*self%dt*aun(ied)*(Vx*dQdx+Vy*dQdy)/(denom)
      apos = max(0.,aun(ied))
      aneg = min(0.,aun(ied))
      flux = Q(p1)*apos + Q(p2)*aneg

      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode) * self%dt / self%vol%array(inode,1)
    end do

  end subroutine MPDATA_Solver__mpdata_abs

end module mpdata_module

