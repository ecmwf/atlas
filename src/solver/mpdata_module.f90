
module mpdata_module
  use model_module
  implicit none

  type, public, extends(Solver_class) :: MPDATA_Solver
    class(Field_class), pointer :: vol ! dual mesh volumes
    class(Field_class), pointer :: S ! dual mesh edge-normals
    real(kind=jprb), dimension(:,:), allocatable :: gradQ
    real(kind=jprb) :: dt

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
    real(kind=jprb), intent(in) :: tmax
    integer :: inode
    character(len=1024)        :: filename
    class(Field_class), pointer :: D,Q,R,grad_D,vol,V

    D => self%state%field("depth")
    Q => self%state%field("momentum")
    R => self%state%field("rhs_forcing")
    grad_D => self%state%field("depth_gradient")
    vol => D%function_space%field("dual_volume")
    V => self%state%field("advective_velocity")


    !write (filename, "(A4,I3.3,A2)") "data",self%iter+1,".d"
    !open(73,file=filename,access='sequential',status='unknown')


    self%dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + self%dt
        
    ! Initialise previous time step
    if (self%iter == 0) then
      write(0,*) "First iteration, computing forcing"
      call self%backup_solution()
      call self%compute_advective_velocities()
      call self%compute_Rn()
    end if
    
    call self%add_forcing_to_solution()

    call self%advect_solution()

    !do inode=1,D%size
      !write(73,'(I6, E20.8, E20.8)') inode, V%array(inode,1), V%array(inode,2)
      !write(73,'(I6, E20.8, E20.8)') inode, R%array(inode,1), R%array(inode,2)
    !  write(73,'(I6, E20.8, E20.8, E20.8)') inode, D%array(inode,1), Q%array(inode,1), Q%array(inode,2)
    !end do
    !close(73)
    

    call self%implicit_solve()

    call self%compute_advective_velocities()
    call self%backup_solution()

    
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

  subroutine MPDATA_Solver__compute_gradient(self,Q,gradQ,Q_is_vector)
    class(MPDATA_Solver), intent(inout) :: self
    real(kind=jprb), dimension(:),   intent(in)    :: Q
    real(kind=jprb), dimension(:,:), intent(inout) :: gradQ
    logical, intent(in) :: Q_is_vector
    
    real(kind=jprb)    :: sx,sy,avgQ
    integer :: ip,iface,p1,p2, f

    ! derivatives 
    gradQ(:,:) = 0
    do iface = 1,self%model%grid%nb_internal_faces
      f = self%model%grid%internal_faces(iface)
      p1  = self%model%grid%faces(f,1)
      p2  = self%model%grid%faces(f,2)
      sx  = self%S%array(f,1)
      sy  = self%S%array(f,2)
      avgQ = ( Q(p1) + Q(p2) )*0.5_jprb
      gradQ(p1,1) = gradQ(p1,1) + sx*avgQ
      gradQ(p2,1) = gradQ(p2,1) - sx*avgQ
      gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
      gradQ(p2,2) = gradQ(p2,2) - sy*avgQ
    end do

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    if (.not. Q_is_vector) then
      do iface = 1,self%model%grid%nb_pole_faces
        f = self%model%grid%pole_faces(iface)
        p1  = self%model%grid%faces(f,1)
        p2  = self%model%grid%faces(f,2)
        sy  = self%S%array(f,2)
        avgQ = ( Q(p1) + Q(p2) )*0.5_jprb
        gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
        gradQ(p2,2) = gradQ(p2,2) + sy*avgQ ! not minus because point at other side of pole
      end do
    end if

  end subroutine MPDATA_Solver__compute_gradient

  subroutine MPDATA_Solver__compute_weighted_average(self,Q,avgQ,absS)
    class(MPDATA_Solver), intent(inout) :: self
    real(kind=jprb), dimension(:),   intent(in)    :: Q
    real(kind=jprb), dimension(:),   intent(inout) :: avgQ
    real(kind=jprb), dimension(:),   intent(inout) :: absS

    real(kind=jprb)    :: sx,sy,Qf, normS
    integer :: ip,iface, f,p1,p2, nb_internal_edges, nb_edges

    nb_internal_edges = self%model%grid%nb_faces
    nb_edges = self%model%grid%nb_faces

    ! derivatives 
    avgQ = 0
    absS = 0

    do iface = 1,self%model%grid%nb_internal_faces
      f = self%model%grid%internal_faces(iface)
      p1 = self%model%grid%faces(f,1)
      p2 = self%model%grid%faces(f,2)
      sx  = self%S%array(f,1)
      sy  = self%S%array(f,2)
      Qf = 0.5*( Q(p1) + Q(p2) )
      normS = sqrt(sx*sx+sy*sy)
      avgQ(p1) = avgQ(p1) + normS*Qf
      avgQ(p2) = avgQ(p2) + normS*Qf
      absS(p1) = absS(p1) + normS
      absS(p2) = absS(p2) + normS
    end do
  end subroutine MPDATA_Solver__compute_weighted_average


  subroutine MPDATA_Solver__mpdata_gauge(self,Q,V,Q_is_vector)
    class(MPDATA_Solver), intent(inout)  :: self
    real(kind=jprb), dimension(:),   intent(inout)  :: Q
    real(kind=jprb), dimension(:,:), intent(in)  :: V
    logical, intent(in) :: Q_is_vector
    
    integer :: inode, iface, pass, p1,p2, nb_nodes, f
    real(kind=jprb) :: advection(self%model%grid%nb_nodes)
    real(kind=jprb) :: aun(self%model%grid%nb_faces)

    real(kind=jprb) :: Qmax(self%model%grid%nb_nodes)
    real(kind=jprb) :: Qmin(self%model%grid%nb_nodes)
    real(kind=jprb) :: rhin(self%model%grid%nb_nodes)
    real(kind=jprb) :: rhout(self%model%grid%nb_nodes)
    real(kind=jprb) :: cp(self%model%grid%nb_nodes)
    real(kind=jprb) :: cn(self%model%grid%nb_nodes)

    real(kind=jprb) :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, eps, apos, aneg

    logical :: limiter = .True.

    eps = 1e-10
    nb_nodes = self%model%grid%nb_nodes

    ! 1. First pass
    ! -------------
      
    ! Compute the rhs in vertices, and normal velocity in faces
    advection = 0.
    do iface = 1,self%model%grid%nb_internal_faces
      f  = self%model%grid%internal_faces(iface)
      p1 = self%model%grid%faces(f,1)
      p2 = self%model%grid%faces(f,2)
      sx = self%S%array(f,1)
      sy = self%S%array(f,2)
      Vx = 0.5*(V(p1,1)+V(p2,1))
      Vy = 0.5*(V(p1,2)+V(p2,2))
      aun(f) = Vx*sx +Vy*sy
      apos = max(0.,aun(f))
      aneg = min(0.,aun(f))
      flux = Q(p1)*apos + Q(p2)*aneg
      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode)/self%vol%array(inode,1) * self%dt
    end do


    ! 2. Second pass
    ! -------------
    
    ! Compute derivatives for mpdata
    call self%compute_gradient(Q,self%gradQ, Q_is_vector )

    ! Compute antidiffusive normal velocity in faces
    do iface = 1,self%model%grid%nb_internal_faces
      f = self%model%grid%internal_faces(iface)
      p1 = self%model%grid%faces(f,1)
      p2 = self%model%grid%faces(f,2)

      ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
      volume_of_two_cells = self%vol%array(p1,1) + self%vol%array(p2,1)
      dQdx = (self%gradQ(p1,1)+self%gradQ(p2,1)) / volume_of_two_cells
      dQdy = (self%gradQ(p1,2)+self%gradQ(p2,2)) / volume_of_two_cells
      Vx = 0.5*(V(p1,1)+V(p2,1))
      Vy = 0.5*(V(p1,2)+V(p2,2))

      ! variable sign option with asymptotic analysis, (mpdata gauge)
      aun(f) = abs(aun(f))*(Q(p2)-Q(p1))*0.5 - 0.5*self%dt*aun(f)*(Vx*dQdx+Vy*dQdy)
    end do

    ! non-isscilatory option (fct)
    if (limiter) then

      Qmax(:)=-1.e10
      Qmin(:)= 1.e10

      do iface = 1,self%model%grid%nb_internal_faces
        f = self%model%grid%internal_faces(iface)
        p1 = self%model%grid%faces(f,1)
        p2 = self%model%grid%faces(f,2)
        Qmax(p1)=max(Qmax(p1),Q(p1),Q(p2))
        Qmin(p1)=min(Qmin(p1),Q(p1),Q(p2))
        Qmax(p2)=max(Qmax(p2),Q(p1),Q(p2))
        Qmin(p2)=min(Qmin(p2),Q(p1),Q(p2))
      enddo
      !write(0,*) "min = ", Qminmin
      !write(0,*) "max = ", Qmaxmax

      do iface = 1,self%model%grid%nb_pole_faces
        f = self%model%grid%pole_faces(iface)
        p1 = self%model%grid%faces(f,1)
        p2 = self%model%grid%faces(f,2)
        if (Q_is_vector) then
          Qmax(p1)=max(Qmax(p1),Q(p1),-Q(p2))
          Qmin(p1)=min(Qmin(p1),Q(p1),-Q(p2))
          Qmax(p2)=max(Qmax(p2),-Q(p1),Q(p2))
          Qmin(p2)=min(Qmin(p2),-Q(p1),Q(p2))
        else
          Qmax(p1)=max(Qmax(p1),Q(p1),Q(p2))
          Qmin(p1)=min(Qmin(p1),Q(p1),Q(p2))
          Qmax(p2)=max(Qmax(p2),Q(p1),Q(p2))
          Qmin(p2)=min(Qmin(p2),Q(p1),Q(p2))
        end if
      enddo

      rhin(:)=0.
      rhout(:)=0.
      cp(:)=0.
      cn(:)=0.

      do iface = 1,self%model%grid%nb_internal_faces
        f = self%model%grid%internal_faces(iface)
        p1 = self%model%grid%faces(f,1)
        p2 = self%model%grid%faces(f,2)
        apos = max(0.,aun(f))
        aneg = min(0.,aun(f))
        rhin(p1)  = rhin(p1) -aneg !+-min(0.,aun(f))
        rhout(p1) = rhout(p1)+ apos !max(0.,aun(f))
        rhin(p2)  = rhin(p2)+ apos !max(0.,aun(f))
        rhout(p2) = rhout(p2)-aneg
      end do

      do inode=1,self%model%grid%nb_nodes
        cp(inode) = ( Qmax(inode)-Q(inode) )*self%vol%array(inode,1)/( rhin(inode) * self%dt + eps )
        cn(inode) = ( Q(inode)-Qmin(inode) )*self%vol%array(inode,1)/( rhout(inode)* self%dt + eps )
      enddo

     ! limited antidiffusive  velocities:
      do iface = 1,self%model%grid%nb_internal_faces
        f = self%model%grid%internal_faces(iface)
        p1 = self%model%grid%faces(f,1)
        p2 = self%model%grid%faces(f,2)
        if(aun(f) > 0.) then
          aun(f)=aun(f)*min(1.,cp(p2),cn(p1))
        else
          aun(f)=aun(f)*min(1.,cn(p2),cp(p1))
        endif
      enddo
    endif

    advection(:) = 0.
    do iface = 1,self%model%grid%nb_internal_faces
      f = self%model%grid%internal_faces(iface)
      p1 = self%model%grid%faces(f,1)
      p2 = self%model%grid%faces(f,2)
      flux = aun(f)
      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode)/self%vol%array(inode,1) * self%dt
    end do

  end subroutine MPDATA_Solver__mpdata_gauge


  subroutine MPDATA_Solver__mpdata_abs(self,Q,V, Q_is_vector)
    class(MPDATA_Solver), intent(inout)  :: self
    real(kind=jprb), dimension(:),   intent(inout)  :: Q
    real(kind=jprb), dimension(:,:), intent(in)  :: V
    logical, intent(in) :: Q_is_vector
    
    integer :: inode, f, iface, pass, p1,p2, apos, aneg, nb_nodes
    real(kind=jprb) :: advection(self%model%grid%nb_nodes)
    real(kind=jprb) :: aun(self%model%grid%nb_faces)
    real(kind=jprb) :: Qabs(self%model%grid%nb_nodes)
    real(kind=jprb) :: Sabs(self%model%grid%nb_nodes)
    real(kind=jprb) :: Qavg(self%model%grid%nb_nodes)

    real(kind=jprb) :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, eps
    real(kind=jprb) :: denom

    eps = 1e-10
    nb_nodes = self%model%grid%nb_nodes

    ! 1. First pass
    ! -------------
      
    ! Compute the rhs in vertices, and normal velocity in faces
    advection = 0.
    do iface = 1,self%model%grid%nb_internal_faces
      f = self%model%grid%internal_faces(iface)
      p1 = self%model%grid%faces(f,1)
      p2 = self%model%grid%faces(f,2)
      sx = self%S%array(f,1)
      sy = self%S%array(f,2)
      Vx = 0.5*(V(p1,1)+V(p2,1))
      Vy = 0.5*(V(p1,2)+V(p2,2))
      aun(f) = Vx*sx + Vy*sy
      apos = max(0.,aun(f))
      aneg = min(0.,aun(f))
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
    call self%compute_gradient(Qabs,self%gradQ, Q_is_vector )
    call self%compute_weighted_average(Q,Qavg,Sabs)


    ! Compute antidiffusive normal velocity in faces
    advection = 0.
    do iface = 1,self%model%grid%nb_internal_faces
      f = self%model%grid%internal_faces(iface)
      p1 = self%model%grid%faces(f,1)
      p2 = self%model%grid%faces(f,2)
      volume_of_two_cells = self%vol%array(p1,1) + self%vol%array(p2,1)
      dQdx = (self%gradQ(p1,1)+self%gradQ(p2,1)) / volume_of_two_cells
      dQdy = (self%gradQ(p1,2)+self%gradQ(p2,2)) / volume_of_two_cells
      Vx = 0.5*(V(p1,1)+V(p2,1))
      Vy = 0.5*(V(p1,2)+V(p2,2))

      ! variable sign option with absolute values
      denom = ( Qavg(p1)+Qavg(p2) ) / ( Sabs(p1) + Sabs(p2) )
      aun(f) = abs(aun(f))*( Qabs(p2)-Qabs(p1) )/( Qabs(p2)+Qabs(p1)+eps ) - 0.5*self%dt*aun(f)*(Vx*dQdx+Vy*dQdy)/(denom)
      apos = max(0.,aun(f))
      aneg = min(0.,aun(f))
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

