!#define OLD_WAY

!#define NDEBUG
#ifndef OLD_WAY
#define NEW_WAY
#endif

#ifndef NDEBUG
#define assert(ASSERTION) \
if( .not. (ASSERTION) ) then; \
write(0,*) "Assertion failed: ASSERTION"; \
call abort; \
endif;
#else
#define assert(ASSERTION)
#endif

#ifndef NDEBUG
#define assert_msg(ASSERTION,msg) \
if( .not. (ASSERTION) ) then; \
write(0,*) "Assertion failed: "//trim(msg); \
call abort; \
endif;
#else
#define assert_msg(ASSERTION,msg)
#endif

#ifdef OLD_WAY
#define LOOP_EDGES \
jedge = 1,geom%nb_internal_edges; \
iedge = geom%internal_edges(jedge);
#define VELOCITY_X  (V(ip1,XX)+V(ip2,XX))*0.5_jprb
#define VELOCITY_Y  (V(ip1,YY)+V(ip2,YY))*0.5_jprb
#else
#define LOOP_EDGES \
iedge = 1,geom%nb_edges;
#define VELOCITY_X  V(iedge,XX)
#define VELOCITY_Y  V(iedge,YY)
#endif


! ===================================================================
! mpdata_module
! -------------
! This module contains strictly algorithmic subroutines
! - mpdata_gauge       : Non-oscillatory variable sign MPDATA advection
! - compute_gradient   : Compute gradient of a scalar array
! ===================================================================
module mpdata_module
  use common_module
  use datastruct_module

  implicit none
  private
  public :: mpdata_gauge
  public :: compute_gradient

  real(kind=jprb), parameter :: eps  = 1.e-10

contains

  subroutine mpdata_gauge(dt,Q,V,order,limited,Q_is_vector,geom)
    real(kind=jprb), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), intent(inout) :: Q(:)
    real(kind=jprb), intent(in) :: V(:,:)
    logical, intent(in) :: Q_is_vector, limited
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, ip1,ip2
    real(kind=jprb) :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, apos, aneg, Q1, Q2
    real(kind=jprb) :: Qmin(geom%nb_nodes)
    real(kind=jprb) :: Qmax(geom%nb_nodes)
    real(kind=jprb) :: rhin(geom%nb_nodes)
    real(kind=jprb) :: rhout(geom%nb_nodes)
    real(kind=jprb) :: cp(geom%nb_nodes)
    real(kind=jprb) :: cn(geom%nb_nodes)
    real(kind=jprb) :: adv(geom%nb_nodes)
    real(kind=jprb) :: aun(geom%nb_edges)
    real(kind=jprb) :: gradQ(geom%nb_nodes,2)
    real(kind=jprb) :: tmp_min, tmp_max

    associate( nb_nodes          => geom%nb_nodes, &
      &        nb_edges          => geom%nb_edges, &
      &        nb_internal_edges => geom%nb_internal_edges, &
      &        nb_pole_edges     => geom%nb_pole_edges, &
      &        internal_edges    => geom%internal_edges, &
      &        pole_edges        => geom%pole_edges, &
      &        vol               => geom%dual_volumes, &
      &        S                 => geom%dual_normals, &
      &        edges             => geom%edges )

    ! 1. First pass
    ! -------------
    
    ! Compute the normal velocity in faces, and advection in vertices
    adv(:) = 0.

    do LOOP_EDGES
      ip1 = edges(iedge,1)
      ip2 = edges(iedge,2)
      Sx = S(iedge,XX)
      Sy = S(iedge,YY)
      Vx = VELOCITY_X
      Vy = VELOCITY_Y
      aun(iedge) = Vx*Sx + Vy*Sy
      apos = max(0._jprb,aun(iedge))
      aneg = min(0._jprb,aun(iedge))
      flux = Q(ip1)*apos + Q(ip2)*aneg
      adv(ip1) = adv(ip1) + flux
      adv(ip2) = adv(ip2) - flux
    end do

    ! Update the unknowns in vertices
    do jnode=1,nb_nodes
      Q(jnode) = Q(jnode) - adv(jnode)/vol(jnode) * dt
    end do

    !write(0,*) "advected = ", L2norm(Q)


    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order
      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector, geom)

      ! Compute antidiffusive normal velocity in faces

      do LOOP_EDGES
        ip1 = edges(iedge,1)
        ip2 = edges(iedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(ip1) + vol(ip2)
        dQdx = (gradQ(ip1,XX)+gradQ(ip2,XX)) / volume_of_two_cells
        dQdy = (gradQ(ip1,YY)+gradQ(ip2,YY)) / volume_of_two_cells
        Vx = VELOCITY_X
        Vy = VELOCITY_Y
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(iedge) = abs(aun(iedge))*(Q(ip2)-Q(ip1))*0.5_jprb &
          &      -0.5_jprb*dt*aun(iedge)*(Vx*dQdx+Vy*dQdy)
      end do

      ! non-oscillatory option
      if (limited) then

        Qmax(:)  = -1.e10
        Qmin(:)  =  1.e10
        rhin(:)  =  0.
        rhout(:) =  0.
        adv(:)   =  0.

        do LOOP_EDGES
          ip1 = edges(iedge,1)
          ip2 = edges(iedge,2)
          Q1 = Q(ip1)
          Q2 = Q(ip2)
          Qmax(ip1) = max( Qmax(ip1), Q1, Q2 )
          Qmin(ip1) = min( Qmin(ip1), Q1, Q2 )
          Qmax(ip2) = max( Qmax(ip2), Q1, Q2 )
          Qmin(ip2) = min( Qmin(ip2), Q1, Q2 )
        end do

        if (Q_is_vector) then
          do jedge = 1,nb_pole_edges
            iedge  = pole_edges(jedge)
            ip1 = edges(iedge,1)
            ip2 = edges(iedge,2)
            Q1 = Q(ip1)
            Q2 = Q(ip2)
            Qmax(ip1) = max( Qmax(ip1), Q1, abs(Q2) )
            Qmin(ip1) = min( Qmin(ip1), Q1,-abs(Q2) )
            Qmax(ip2) = max( Qmax(ip2), Q2, abs(Q1) )
            Qmin(ip2) = min( Qmin(ip2), Q2,-abs(Q1) )
          end do
        end if

        do LOOP_EDGES
          ip1 = edges(iedge,1)
          ip2 = edges(iedge,2)
          apos = max(0._jprb,aun(iedge))
          aneg = min(0._jprb,aun(iedge))
          rhin(ip1)  = rhin(ip1)  - aneg
          rhout(ip1) = rhout(ip1) + apos
          rhin(ip2)  = rhin(ip2)  + apos
          rhout(ip2) = rhout(ip2) - aneg
        end do

        do jnode=1,nb_nodes
          assert(vol(jnode) .ne. 0._jprb )
          cp(jnode) = ( Qmax(jnode)-Q(jnode) )*vol(jnode)/( rhin(jnode) * dt + eps )
          cn(jnode) = ( Q(jnode)-Qmin(jnode) )*vol(jnode)/( rhout(jnode)* dt + eps )
        enddo

       ! limited antidiffusive  velocities:

        do LOOP_EDGES
          ip1 = edges(iedge,1)
          ip2 = edges(iedge,2)
          if(aun(iedge) > 0.) then
            aun(iedge)=aun(iedge)*min(1._jprb,cp(ip2),cn(ip1))
          else
            aun(iedge)=aun(iedge)*min(1._jprb,cn(ip2),cp(ip1))
          endif
        enddo
      endif

      ! Compute fluxes from (limited) antidiffusive velocity
      do LOOP_EDGES
        ip1 = edges(iedge,1)
        ip2 = edges(iedge,2)
        flux = aun(iedge)
        adv(ip1) = adv(ip1) + flux
        adv(ip2) = adv(ip2) - flux
      end do

      ! Update the unknowns in vertices
      do jnode=1,nb_nodes
        Q(jnode) = Q(jnode) - adv(jnode)/vol(jnode) * dt
      end do

    end do ! other passes
    end associate
  end subroutine mpdata_gauge

  subroutine compute_gradient(Q,gradQ,Q_is_vector,geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), intent(in)    :: Q(:)
    real(kind=jprb), intent(inout) :: gradQ(:,:)
    logical, intent(in) :: Q_is_vector
    
    real(kind=jprb) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2

    ! derivatives 
    gradQ(:,:) = 0

    do LOOP_EDGES
      ip1  = geom%edges(iedge,1)
      ip2  = geom%edges(iedge,2)
      Sx  = geom%dual_normals(iedge,XX)
      Sy  = geom%dual_normals(iedge,YY)
      avgQ = ( Q(ip1) + Q(ip2) )*0.5_jprb
      gradQ(ip1,XX) = gradQ(ip1,XX) + Sx*avgQ
      gradQ(ip2,XX) = gradQ(ip2,XX) - Sx*avgQ
      gradQ(ip1,YY) = gradQ(ip1,YY) + Sy*avgQ
      gradQ(ip2,YY) = gradQ(ip2,YY) - Sy*avgQ
    end do

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    if (.not. Q_is_vector) then
      do jedge = 1,geom%nb_pole_edges
        iedge    = geom%pole_edges(jedge)
        ip1   = geom%edges(iedge,1)
        ip2   = geom%edges(iedge,2)
        Sy   = geom%dual_normals(iedge,YY)
        avgQ = ( Q(ip1) + Q(ip2) )*0.5_jprb
        ! correct for wrongly subtracting sy*avgQ in previous loop,
        ! instead of adding, because point at other side of pole
#ifdef  OLD_WAY
        gradQ(ip1,YY) = gradQ(ip1,YY) + Sy*avgQ
        gradQ(ip2,YY) = gradQ(ip2,YY) + Sy*avgQ
#else
        gradQ(ip2,YY) = gradQ(ip2,YY) + 2*Sy*avgQ 
#endif
      end do
    end if
  end subroutine compute_gradient

end module mpdata_module


! ===================================================================
! shallow_water_module
! --------------------
! This module contains subroutines to solve the shallow water eqs
! on a sphere
! - set_time_step : To set the solver's time step
! - compute_metrics : To modify dual-volumes to map to the sphere
! - allocations   : To pre-allocate necessary structures and states
! - set_state_rossby_haurwits : To initialise the state with
!                               Rossby Haurwitz waves
! - propagate_state : To propagate the state with a given time step
! ===================================================================
module shallow_water_module
  use common_module
  use datastruct_module
  use mpdata_module, &
    & only: mpdata_gauge, &
    &       compute_gradient

  implicit none
  private
  public :: compute_metrics 
  public :: propagate_state
  public :: allocations
  public :: set_state_rossby_haurwitz
  public :: set_time_step

  real(kind=jprb), parameter :: eps    = 1.e-10
  real(kind=jprb), parameter :: radius = 6371.22e+03
  real(kind=jprb), parameter :: f0     = 1.4584e-04 !coriolis parameter (=2xearth's omega)
  real(kind=jprb), parameter :: grav   = 9.80616
  real(kind=jprb), parameter :: pi     = acos(-1._jprb)


  real(kind=jprb) :: dt_forward = 20.
  integer :: iter = 0

contains
 
  subroutine set_time_step(dt)
    real(kind=jprb), intent(in) :: dt
    dt_forward = dt
  end subroutine set_time_step

  subroutine compute_metrics(geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb) :: y, G, cos_y, sin_y, hx, hy
    integer :: jnode
    do jnode=1,geom%nb_nodes
      y = geom%coordinates(jnode,YY)
      cos_y = cos(y)
      sin_y = sin(y)
      hx = radius*cos_y
      hy = radius
      G = hx*hy
      geom%dual_volumes(jnode) = geom%dual_volumes(jnode)*G
    enddo
  end subroutine compute_metrics



  subroutine allocations(geom)
    type(DataStructure_type), intent(inout) :: geom
    call create_state_fields(geom)
  end subroutine allocations



  subroutine propagate_state(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout), target :: geom
    real(kind=jprb) :: tmax, t0, dt_fwd
    character(len=200) :: step_info
    tmax = geom%fields%time+dt
    do while (geom%fields%time < tmax)
      t0 = geom%fields%time
      dt_fwd = min( dt_forward, tmax-t0 )
      call step_forward(iter,dt_fwd,geom)

      write(step_info,'(A6,I8,A12,F9.1,A12,F8.1)') "step = ",iter, &
         & "  time = ",geom%fields%time, &
         & "  dt = ",dt_fwd
      call log_info( step_info )
    end do
  end subroutine propagate_state




  subroutine step_forward(step,dt,geom)
    integer, intent(inout) :: step
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom

    if (step == 0) then ! Pre-compute forcing
      call backup_solution(geom)
      call compute_advective_velocities(geom)
      call compute_forcing(geom)
    end if
    
    call add_forcing_to_solution(dt,geom)
    call advect_solution(dt,geom)
    call implicit_solve(dt,geom)
    call compute_advective_velocities(geom)
    call backup_solution(geom)

    geom%fields%time = geom%fields%time+dt
    step = step+1

  end subroutine step_forward




  subroutine create_state_fields(geom)
    type(DataStructure_type), intent(inout) :: geom
    call create_scalar_field_in_nodes("depth",geom)
    call create_vector_field_in_nodes("momentum",geom)
    call create_vector_field_in_nodes("momentum_forcing",geom)
    call create_scalar_field_in_nodes("depth_backup",geom)
    call create_vector_field_in_nodes("momentum_backup",geom)
#ifdef OLD_WAY
    call create_vector_field_in_nodes("advective_velocity",geom)
#else
    call create_vector_field_in_edges("advective_velocity",geom)
#endif
  end subroutine create_state_fields



  subroutine set_state_rossby_haurwitz(geom)
    type(DataStructure_type), intent(inout)      :: geom
    real(kind=jprb), dimension(:), pointer   :: D
    real(kind=jprb), dimension(:,:), pointer :: Q
    integer :: jnode
    integer :: ir,ip
    real(kind=jprb)    :: aaa0,zk,om,ph0,g,ath,bth,cth,x,y,th,cor

    ! statement-functions, before first statement
    ATH(TH) = om*0.5*(f0+om)*(cos(TH))**2 &
      & +0.25*zk**2*(cos(TH))**(2*ir)*( (ir+1)*(cos(TH))**2 &
      & +float(2*ir**2-ir-2)-2.*ir**2/(cos(TH))**2 )
    BTH(TH) = (f0+2.*om)*zk/float((ir+1)*(ir+2))*(cos(TH))**ir &
      & *( float(ir**2+2*ir+2)-((ir+1)*cos(TH))**2 )
    CTH(TH) = 0.25*zk**2*(cos(TH))**(2*ir)*( float(ir+1)*(cos(TH))**2 &
      & -float(ir+2) )  

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3
    aaa0 = 0.

    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)

    do jnode=1,geom%nb_nodes
      x=geom%coordinates(jnode,XX)
      y=geom%coordinates(jnode,YY)
      if(x == 2._jprb*pi) x=0.
      cor=f0*sin(y)
      Q(jnode,XX) =  radius*OM*cos(y)+radius*ZK*cos(IR*x) *(cos(y))**(IR-1)*(IR*(sin(y))**2-(cos(y))**2)
      Q(jnode,YY) = -radius*ZK*IR*(cos(y))**(IR-1)*sin(y)*sin(IR*x)
      D(jnode) = (ph0+radius**2*ATH(y)+radius**2*BTH(y)*cos(IR*x)+radius**2*CTH(y)*cos(2._jprb*IR*x)) &
          & /(grav)
      D(jnode) = max(aaa0,D(jnode))
      Q(jnode,XX) = Q(jnode,XX) * D(jnode)
      Q(jnode,YY) = Q(jnode,YY) * D(jnode)
      if(y == 0.5_jprb*pi) Q(jnode,XX)=0.
      if(y ==-0.5_jprb*pi) Q(jnode,XX)=0.
    end do
  end subroutine set_state_rossby_haurwitz



  subroutine backup_solution(geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), dimension(:),   pointer :: D, D0
    real(kind=jprb), dimension(:,:), pointer :: Q, Q0
    D  => scalar_field("depth",geom)
    D0 => scalar_field("depth_backup",geom)
    Q  => vector_field("momentum",geom)
    Q0 => vector_field("momentum_backup",geom)

    D0(:)   = D(:)
    Q0(:,:) = Q(:,:)
  end subroutine backup_solution



  subroutine compute_advective_velocities(geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb) :: y, y1,y2, hx, hy, Qx, Qy, Q0x, Q0y, Dmod, D0mod
    integer :: jnode, jedge, iedge, ip1, ip2
    real(kind=jprb), dimension(:),   pointer :: D, D0
    real(kind=jprb), dimension(:,:), pointer :: Q, Q0
#ifdef OLD_WAY
    real(kind=jprb), dimension(:,:), pointer :: Vnodes
    Vnodes  => vector_field("advective_velocity",geom)
#else
    real(kind=jprb), dimension(:,:), pointer :: Vedges
    real(kind=jprb)                          :: Vnodes(geom%nb_nodes,2)
    Vedges  => vector_field("advective_velocity",geom)
#endif
    D  => scalar_field("depth",geom)
    D0 => scalar_field("depth_backup",geom)
    Q  => vector_field("momentum",geom)
    Q0 => vector_field("momentum_backup",geom)
    do jnode=1,geom%nb_nodes
      Qx    = Q(jnode,XX)
      Qy    = Q(jnode,YY)
      Dmod  = max( eps, D(jnode) )
      Q0x   = Q0(jnode,XX)
      Q0y   = Q0(jnode,YY)
      D0mod = max( eps, D0(jnode) )
      
      y  = geom%coordinates(jnode,YY)
      hx = radius*cos(y)
      hy = radius
      ! this really computes V = G*contravariant_velocity, 
      ! with    G=hx*hy,
      !         physical_velocity = dotproduct( [hx,hy] , contravariant_velocity )
      ! V = (hx*hy) * [u/hx, v/hy] = [u*hy, v*hx]
      ! and hx = r*cos(y)  ,  hy = r
      ! and Q = [ D*u , D*v ]
      Vnodes(jnode,XX) = ( 1.5_jprb*Qx/Dmod - 0.5_jprb*Q0x/D0mod ) * hy
      Vnodes(jnode,YY) = ( 1.5_jprb*Qy/Dmod - 0.5_jprb*Q0y/D0mod ) * hx
    end do

#ifdef NEW_WAY
    do jedge=1,geom%nb_edges
      ip1 = geom%edges(jedge,1)
      ip2 = geom%edges(jedge,2)
      ! Evaluate hx and hy at the edge midpoint
      y1 = geom%coordinates(ip1,YY)
      y2 = geom%coordinates(ip2,YY)
      y = (y1+y2)*0.5_jprb
      hx    = radius*cos(y)
      hy    = radius
      Vedges(jedge,XX) = (Vnodes(ip1,XX)+Vnodes(ip2,XX))*0.5_jprb
      Vedges(jedge,YY) = (Vnodes(ip1,YY)+Vnodes(ip2,YY))*0.5_jprb
    enddo
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,geom%nb_pole_edges
      iedge = geom%pole_edges(jedge)
      Vedges(iedge,YY) = 0.
    enddo
#endif
  end subroutine compute_advective_velocities

  subroutine compute_forcing(geom)
    type(DataStructure_type), intent(inout) :: geom
    integer :: jnode
    real(kind=jprb) :: cor, y, sin_y, cos_y, Qx, Qy, Dmod, dDdx, dDdy, vol, hx, hy
    real(kind=jprb), dimension(:),   pointer :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    real(kind=jprb) :: grad_D(geom%nb_nodes, 2)
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)

    call compute_gradient( D, grad_D, .False., geom )
    
    do jnode=1,geom%nb_nodes
      y     = geom%coordinates(jnode,YY)
      sin_y = sin(y)
      cos_y = max( eps, cos(y) )
      Qx    = Q(jnode,XX)
      Qy    = Q(jnode,YY)
      Dmod  = max( eps, D(jnode) )
      vol   = geom%dual_volumes(jnode)
      hx    = radius*cos_y
      hy    = radius
      cor   = f0 * sin_y
      R(jnode,XX) = -grav*D(jnode)*grad_D(jnode,XX)*hy/vol + cor*Qy + sin_y/(radius*cos_y)*Qx*Qy/Dmod
      R(jnode,YY) = -grav*D(jnode)*grad_D(jnode,YY)*hx/vol - cor*Qx - sin_y/(radius*cos_y)*Qx*Qx/Dmod
    end do

  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    integer :: jnode, nb_nodes
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)
    do jnode=1,geom%nb_nodes
      Q(jnode,XX) = Q(jnode,XX) + 0.5_jprb*dt*R(jnode,XX)
      Q(jnode,YY) = Q(jnode,YY) + 0.5_jprb*dt*R(jnode,YY)
    end do
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    integer :: jnode, m
    real(kind=jprb) :: cor, y, sin_y, cos_y, Qx, Qy, Rx, Ry, Dmod, dDdx, dDdy, hx, hy
    real(kind=jprb) :: Qx_adv, Qy_adv, Rx_exp, Ry_exp, vol

    real(kind=jprb), dimension(:),   pointer :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    real(kind=jprb) :: grad_D(geom%nb_nodes, 2)

    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    call compute_gradient( D, grad_D, .False., geom )
    
    do jnode=1,geom%nb_nodes
      y     = geom%coordinates(jnode,YY)
      sin_y = sin(y)
      cos_y = max( eps, cos(y) )
      Qx    = Q(jnode,XX)
      Qy    = Q(jnode,YY)
      Dmod  = max( eps, D(jnode) )
      vol   = geom%dual_volumes(jnode)
      hx    = radius*cos_y
      hy    = radius
      cor   = f0 * sin_y

      Rx_exp = -grav*D(jnode)*grad_D(jnode,XX)*hy/vol
      Ry_exp = -grav*D(jnode)*grad_D(jnode,YY)*hx/vol

      Qx_adv = Qx
      Qy_adv = Qy

      do m=1,3 ! Three iterations at most is enough to converge
        Rx = Rx_exp + cor*Qy + sin_y/(radius*cos_y)*Qx*Qy/Dmod
        Ry = Ry_exp - cor*Qx - sin_y/(radius*cos_y)*Qx*Qx/Dmod
        Qx = Qx_adv + 0.5_jprb*dt*Rx
        Qy = Qy_adv + 0.5_jprb*dt*Ry
      end do
      Q(jnode,XX) = Qx
      Q(jnode,YY) = Qy
      R(jnode,XX) = Rx_exp + cor*Qy + sin_y/(radius*cos_y)*Qx*Qy/Dmod
      R(jnode,YY) = Ry_exp - cor*Qx - sin_y/(radius*cos_y)*Qx*Qx/Dmod
    end do
  end subroutine implicit_solve



  subroutine advect_solution(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), dimension(:),   pointer :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, V
    integer :: jnode, ip1, ip2
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    V => vector_field("advective_velocity",geom)

    !    mpdata_gauge( time, variable, velocity, order, limit,  is_vector, geom )
    call mpdata_gauge( dt,   D,        V,        2,     .True., .False.,   geom )
    call mpdata_gauge( dt,   Q(:,XX),  V,        2,     .True., .True. ,   geom )
    call mpdata_gauge( dt,   Q(:,YY),  V,        2,     .True., .True. ,   geom )

    ! remove noise from periodic boundary.. why is it even there in the mesh?
    do jnode=1,geom%nb_ghost_nodes
      ip1 = geom%ghost_nodes(jnode,1)
      ip2 = geom%ghost_nodes(jnode,2)
      D(ip2)    = D(ip1)
      Q(ip2,XX) = Q(ip1,XX)
      Q(ip2,YY) = Q(ip1,YY)
    end do
  end subroutine advect_solution

end module shallow_water_module



! ===================================================================
! shallow_water program
! ---------------------
! This program solves the shallow water equations on a sphere,
! initialised with the Rossby-Haurwitz waves.
! ===================================================================
program shallow_water

  use common_module
  use gmsh_module, only: write_gmsh_mesh, write_gmsh_state
  use grib_module, only: write_grib

  use read_joanna_module, only: read_joanna
  use datastruct_module,  only: DataStructure_type
  use shallow_water_module, only: &
    & compute_metrics, &
    & allocations, &
    & set_state_rossby_haurwitz, &
    & set_time_step, &
    & propagate_state
  implicit none

  ! Configuration parameters
  real(kind=jprb) :: dt = 20.              ! solver time-step
  integer         :: nb_steps = 15         ! Number of propagations
  integer         :: hours_per_step = 24   !Propagation time
  logical         :: write_itermediate_output = .True.

  ! Declarations
  type(DataStructure_type) :: g
  real(kind=jprb), parameter :: hours = 3600.     ! in seconds
  real(kind=jprb), parameter :: days  = 24.*hours ! in seconds
  integer :: jstep = 0
  character(len=1024) :: filename
  character(len=1024) :: string
  integer :: clck_counts_start, clck_counts_beg, clck_counts_end, clck_rate
 
  call set_log_level(LOG_LEVEL_INFO)
  call log_info("Going to solve shallow_water equations for "//trim(str(nb_steps))//" days")

  ! Execution
  call read_joanna("data/meshvol.d", g)
  call compute_metrics(g)
  call allocations(g)
  call set_state_rossby_haurwitz(g)
  call set_time_step( dt )

  call write_gmsh_mesh(g%internal_mesh,"data/mesh.msh")

  write (filename, "(A,I2.2,A)") "data/fields",jstep,".msh"
  call write_gmsh_state(g%fields,filename)

  write (filename, "(A,I2.2,A)") "data/depth",jstep,".grib"
  call write_grib(g,filename)

  call system_clock ( clck_counts_start, clck_rate )
  do jstep=1,nb_steps 

    call system_clock ( clck_counts_beg, clck_rate )
  
    call propagate_state( hours_per_step*hours, g)
  
    call system_clock ( clck_counts_end, clck_rate )
    write (string, *) "propagated to ",jstep*hours_per_step," hours.   time = ",(clck_counts_end - clck_counts_beg) / real(clck_rate),&
     &  "walltime = ",(clck_counts_end - clck_counts_start) / real(clck_rate), new_line('A')
    call log_info(string)
    
    if (write_itermediate_output) then
      write (filename, "(A,I2.2,A)") "data/fields",jstep,".msh"
      call write_gmsh_state(g%fields,filename)

      write (filename, "(A,I2.2,A)") "data/depth",jstep,".grib"
      call write_grib(g,filename)
    end if
  end do ! steps

  if (.not. write_itermediate_output) then
    write (filename, "(A,I2.2,A)") "data/fields",jstep,".msh"
    call write_gmsh_state(g%fields,filename)

    write (filename, "(A,I2.2,A)") "data/depth",jstep,".grib"
    call write_grib(g,filename)
  end if

end program shallow_water
