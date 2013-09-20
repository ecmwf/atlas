
module mpdata_module
  use model_module
  implicit none

  type, public, extends(Solver_class) :: MPDATA_Solver
    class(Field_class), pointer :: vol ! dual mesh volumes
    class(Field_class), pointer :: S ! dual mesh edge-normals
    real(kind=jprb), dimension(:,:), allocatable :: gradQ
    real(kind=jprb) :: dt

    ! needed for advection
    real(kind=jprb), dimension(:), allocatable, private :: advection
    real(kind=jprb), dimension(:), allocatable, private :: aun
    
    ! needed for limiter
    real(kind=jprb), dimension(:), allocatable, private :: Qmax
    real(kind=jprb), dimension(:), allocatable, private :: Qmin
    real(kind=jprb), dimension(:), allocatable, private :: rhin
    real(kind=jprb), dimension(:), allocatable, private :: rhout
    real(kind=jprb), dimension(:), allocatable, private :: cp
    real(kind=jprb), dimension(:), allocatable, private :: cn

  contains
    procedure, public,  pass :: init => MPDATA_Solver__init
    procedure, public,  pass :: step => MPDATA_Solver__step
    procedure, pass :: compute_gradient => MPDATA_Solver__compute_gradient
    procedure, pass :: compute_forcing   => MPDATA_Solver__compute_forcing
    procedure, pass :: implicit_solve => MPDATA_Solver__implicit_solve
    procedure, pass :: compute_advective_velocities => MPDATA_Solver__compute_advective_velocities
    procedure, pass :: mpdata_gauge => MPDATA_Solver__mpdata_gauge
    procedure, pass :: backup_solution => MPDATA_Solver__backup_solution
    procedure, pass :: add_forcing_to_solution => MPDATA_Solver__add_forcing_to_solution
    procedure, pass :: advect_solution => MPDATA_Solver__advect_solution

    
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

    allocate( self%aun(self%model%grid%nb_faces) )
    allocate( self%advection(self%model%grid%nb_nodes) )
    allocate( self%Qmax(self%model%grid%nb_nodes) )
    allocate( self%Qmin(self%model%grid%nb_nodes) )
    allocate( self%rhin(self%model%grid%nb_nodes) )
    allocate( self%rhout(self%model%grid%nb_nodes) )
    allocate( self%cp(self%model%grid%nb_nodes) )
    allocate( self%cn(self%model%grid%nb_nodes) )

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

    self%dt = min( self%dt_stability, tmax-self%state%time )
    self%state%time = self%state%time + self%dt
        
    ! Initialise previous time step
    if (self%iter == 0) then
      write(0,*) "First iteration, computing forcing"
      call self%backup_solution()
      call self%compute_advective_velocities()
      call self%compute_forcing()
    end if
    
    call self%add_forcing_to_solution()

    call self%advect_solution()

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
    ! Look inside shallow_water_module for specific implementation
  end subroutine MPDATA_Solver__backup_solution

  subroutine MPDATA_Solver__add_forcing_to_solution(self)
    class(MPDATA_Solver), intent(inout) :: self
    ! Look inside shallow_water_module for specific implementation
  end subroutine MPDATA_Solver__add_forcing_to_solution

  subroutine MPDATA_Solver__compute_advective_velocities(self)
    class(MPDATA_Solver), intent(inout)    :: self
    ! Look inside shallow_water_module for specific implementation
  end subroutine MPDATA_Solver__compute_advective_velocities

  subroutine MPDATA_Solver__advect_solution(self)
    class(MPDATA_Solver), intent(inout)    :: self
    ! Look inside shallow_water_module for specific implementation
  end subroutine MPDATA_Solver__advect_solution

  subroutine MPDATA_Solver__compute_forcing(self)
    class(MPDATA_Solver), intent(inout)    :: self
    ! Look inside shallow_water_module for specific implementation
  end subroutine MPDATA_Solver__compute_forcing

  subroutine MPDATA_Solver__implicit_solve(self)
    class(MPDATA_Solver), intent(inout)    :: self
    ! Look inside shallow_water_module for specific implementation
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

    associate( face_nodes        => self%model%grid%faces , &
      &        internal_faces    => self%model%grid%internal_faces , &
      &        pole_faces        => self%model%grid%pole_faces , &
      &        nb_nodes          => self%model%grid%nb_nodes , &
      &        nb_pole_faces     => self%model%grid%nb_pole_faces , &
      &        nb_internal_faces => self%model%grid%nb_internal_faces, &
      &        S                 => self%S%array )

    ! derivatives 
    gradQ(:,:) = 0
    do iface = 1,nb_internal_faces
      f = internal_faces(iface)
      p1  = face_nodes(f,1)
      p2  = face_nodes(f,2)
      sx  = S(f,1)
      sy  = S(f,2)
      avgQ = ( Q(p1) + Q(p2) )*0.5_jprb
      gradQ(p1,1) = gradQ(p1,1) + sx*avgQ
      gradQ(p2,1) = gradQ(p2,1) - sx*avgQ
      gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
      gradQ(p2,2) = gradQ(p2,2) - sy*avgQ
    end do

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    if (.not. Q_is_vector) then
      do iface = 1,nb_pole_faces
        f = pole_faces(iface)
        p1  = face_nodes(f,1)
        p2  = face_nodes(f,2)
        sy  = S(f,2)
        avgQ = ( Q(p1) + Q(p2) )*0.5_jprb
        gradQ(p1,2) = gradQ(p1,2) + sy*avgQ
        gradQ(p2,2) = gradQ(p2,2) + sy*avgQ ! not minus because point at other side of pole
      end do
    end if
    end associate
  end subroutine MPDATA_Solver__compute_gradient


  subroutine MPDATA_Solver__mpdata_gauge(self,Q,V,order,limited,Q_is_vector)
    class(MPDATA_Solver), intent(inout)  :: self
    real(kind=jprb), dimension(:),   intent(inout)  :: Q
    real(kind=jprb), dimension(:,:), intent(in)  :: V
    logical, intent(in) :: Q_is_vector, limited
    integer, intent(in) :: order
    
    integer :: inode, iface, pass, p1,p2, f
    real(kind=jprb) :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, apos, aneg
    real(kind=jprb), parameter :: eps = 1e-10

    associate( face_nodes        => self%model%grid%faces , &
      &        internal_faces    => self%model%grid%internal_faces , &
      &        pole_faces        => self%model%grid%pole_faces , &
      &        nb_nodes          => self%model%grid%nb_nodes , &
      &        nb_pole_faces     => self%model%grid%nb_pole_faces , &
      &        nb_internal_faces => self%model%grid%nb_internal_faces, &
      &        Qmin              => self%Qmin, &
      &        Qmax              => self%Qmax, &
      &        advection         => self%advection, &
      &        aun               => self%aun, &
      &        rhin              => self%rhin, &
      &        rhout             => self%rhout, &
      &        cp                => self%cp, &
      &        cn                => self%cn, &
      &        dt                => self%dt, &
      &        vol               => self%vol%array(:,1),&
      &        S                 => self%S%array, &
      &        gradQ             => self%gradQ )

    ! 1. First pass
    ! -------------
      
    ! Compute the normal velocity in faces, and advection in vertices
    advection(:) = 0.
    do iface = 1,nb_internal_faces
      f  = internal_faces(iface)
      p1 = face_nodes(f,1)
      p2 = face_nodes(f,2)
      sx = S(f,1)
      sy = S(f,2)
      Vx = 0.5_jprb*(V(p1,1)+V(p2,1))
      Vy = 0.5_jprb*(V(p1,2)+V(p2,2))
      aun(f) = Vx*sx +Vy*sy
      apos = max(0._jprb,aun(f))
      aneg = min(0._jprb,aun(f))
      flux = Q(p1)*apos + Q(p2)*aneg
      advection(p1) = advection(p1) + flux
      advection(p2) = advection(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - advection(inode)/vol(inode) * dt
    end do


    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do pass=2,order
      ! Compute derivatives for mpdata
      call self%compute_gradient(Q,gradQ, Q_is_vector )

      ! Compute antidiffusive normal velocity in faces
      do iface = 1,nb_internal_faces
        f  = internal_faces(iface)
        p1 = face_nodes(f,1)
        p2 = face_nodes(f,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(p1) + vol(p2)
        dQdx = (gradQ(p1,1)+gradQ(p2,1)) / volume_of_two_cells
        dQdy = (gradQ(p1,2)+gradQ(p2,2)) / volume_of_two_cells
        Vx = 0.5_jprb*(V(p1,1)+V(p2,1))
        Vy = 0.5_jprb*(V(p1,2)+V(p2,2))

        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(f) = abs(aun(f))*(Q(p2)-Q(p1))*0.5_jprb - 0.5_jprb*dt*aun(f)*(Vx*dQdx+Vy*dQdy)
      end do

      ! non-isscilatory option
      if (limited) then

        Qmax(:)=-1.e10
        Qmin(:)= 1.e10

        do iface = 1,nb_internal_faces
          f  = internal_faces(iface)
          p1 = face_nodes(f,1)
          p2 = face_nodes(f,2)
          Qmax(p1)=max(Qmax(p1),Q(p1),Q(p2))
          Qmin(p1)=min(Qmin(p1),Q(p1),Q(p2))
          Qmax(p2)=max(Qmax(p2),Q(p1),Q(p2))
          Qmin(p2)=min(Qmin(p2),Q(p1),Q(p2))
        enddo

        do iface = 1,nb_pole_faces
          f  = pole_faces(iface)
          p1 = face_nodes(f,1)
          p2 = face_nodes(f,2)
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

        do iface = 1,nb_internal_faces
          f  = internal_faces(iface)
          p1 = face_nodes(f,1)
          p2 = face_nodes(f,2)
          apos = max(0._jprb,aun(f))
          aneg = min(0._jprb,aun(f))
          rhin(p1)  = rhin(p1)  - aneg
          rhout(p1) = rhout(p1) + apos
          rhin(p2)  = rhin(p2)  + apos
          rhout(p2) = rhout(p2) - aneg
        end do

        do inode=1,nb_nodes
          cp(inode) = ( Qmax(inode)-Q(inode) )*vol(inode)/( rhin(inode) * dt + eps )
          cn(inode) = ( Q(inode)-Qmin(inode) )*vol(inode)/( rhout(inode)* dt + eps )
        enddo

       ! limited antidiffusive  velocities:
        do iface = 1,nb_internal_faces
          f  = internal_faces(iface)
          p1 = face_nodes(f,1)
          p2 = face_nodes(f,2)
          if(aun(f) > 0.) then
            aun(f)=aun(f)*min(1._jprb,cp(p2),cn(p1))
          else
            aun(f)=aun(f)*min(1._jprb,cn(p2),cp(p1))
          endif
        enddo
      endif

      ! Compute fluxes from (limited) antidiffusive velocity
      advection(:) = 0.
      do iface = 1,nb_internal_faces
        f  = internal_faces(iface)
        p1 = face_nodes(f,1)
        p2 = face_nodes(f,2)
        flux = aun(f)
        advection(p1) = advection(p1) + flux
        advection(p2) = advection(p2) - flux
      end do

      ! Update the unknowns in vertices
      do inode=1,nb_nodes
        Q(inode) = Q(inode) - advection(inode)/vol(inode) * dt
      end do

    end do ! other passes

    end associate

  end subroutine MPDATA_Solver__mpdata_gauge


end module mpdata_module

