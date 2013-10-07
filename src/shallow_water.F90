#define OLD_WAY
#define NDEBUG
#ifndef NDEBUG
#define assert(ASSERTION) \
if( .not. (ASSERTION) ) then; \
write(0,*) "Assertion failed: ASSERTION"; \
call abort; \
endif;
#else
#define assert(ASSERTION)
#endif

#ifdef OLD_WAY
#define LOOP_EDGES(EDGE) \
iedge = 1,geom%nb_internal_edges; \
EDGE = geom%internal_edges(iedge);
#define VELOCITY_X  (V(p1,XX)+V(p2,XX))*0.5_jprb
#define VELOCITY_Y  (V(p1,YY)+V(p2,YY))*0.5_jprb
#else
#define LOOP_EDGES(EDGE) \
EDGE = 1,geom%nb_edges;
#define VELOCITY_X  V(e,XX)
#define VELOCITY_Y  V(e,YY)
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
    type(Geometry_class), intent(inout) :: geom
    real(kind=jprb), intent(inout) :: Q(:)
    real(kind=jprb), intent(in) :: V(:,:)
    logical, intent(in) :: Q_is_vector, limited
    integer, intent(in) :: order
    integer :: inode, iedge, pass, p1,p2, e
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

    do LOOP_EDGES(e)
      p1 = edges(e,1)
      p2 = edges(e,2)
      Sx = S(e,XX)
      Sy = S(e,YY)
      Vx = VELOCITY_X
      Vy = VELOCITY_Y
      aun(e) = Vx*Sx + Vy*Sy
      apos = max(0._jprb,aun(e))
      aneg = min(0._jprb,aun(e))
      flux = Q(p1)*apos + Q(p2)*aneg
      adv(p1) = adv(p1) + flux
      adv(p2) = adv(p2) - flux
    end do

    ! Update the unknowns in vertices
    do inode=1,nb_nodes
      Q(inode) = Q(inode) - adv(inode)/vol(inode) * dt
    end do


    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do pass=2,1!order
      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector, geom)

      ! Compute antidiffusive normal velocity in faces

      do LOOP_EDGES(e)
        p1 = edges(e,1)
        p2 = edges(e,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(p1) + vol(p2)
        dQdx = (gradQ(p1,XX)+gradQ(p2,XX)) / volume_of_two_cells
        dQdy = (gradQ(p1,YY)+gradQ(p2,YY)) / volume_of_two_cells
        Vx = VELOCITY_X
        Vy = VELOCITY_Y
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(e) = abs(aun(e))*(Q(p2)-Q(p1))*0.5_jprb &
          &      -0.5_jprb*dt*aun(e)*(Vx*dQdx+Vy*dQdy)
      end do

      ! non-oscillatory option
      if (limited) then

        Qmax(:)  = -1.e10
        Qmin(:)  =  1.e10
        rhin(:)  =  0.
        rhout(:) =  0.
        adv(:)   =  0.

        do iedge = 1,nb_internal_edges
          e  = internal_edges(iedge)
          p1 = edges(e,1)
          p2 = edges(e,2)
          Q1 = Q(p1)
          Q2 = Q(p2)
          Qmax(p1)=max(Qmax(p1),Q1,Q2)
          Qmin(p1)=min(Qmin(p1),Q1,Q2)
          Qmax(p2)=max(Qmax(p2),Q1,Q2)
          Qmin(p2)=min(Qmin(p2),Q1,Q2)
        enddo

        if (Q_is_vector) then
          do iedge = 1,nb_pole_edges
            e  = pole_edges(iedge)
            p1 = edges(e,1)
            p2 = edges(e,2)
            Q1 = Q(p1)
            Q2 = Q(p2)
            if (Q_is_vector) then
              Qmax(p1)=max(Qmax(p1),Q1,-Q2)
              Qmin(p1)=min(Qmin(p1),Q1,-Q2)
              Qmax(p2)=max(Qmax(p2),-Q1,Q2)
              Qmin(p2)=min(Qmin(p2),-Q1,Q2)
            else
              Qmax(p1)=max(Qmax(p1),Q1,Q2)
              Qmin(p1)=min(Qmin(p1),Q1,Q2)
              Qmax(p2)=max(Qmax(p2),Q1,Q2)
              Qmin(p2)=min(Qmin(p2),Q1,Q2)
            endif
          enddo
        end if
     

      do LOOP_EDGES(e)
          p1 = edges(e,1)
          p2 = edges(e,2)
          apos = max(0._jprb,aun(e))
          aneg = min(0._jprb,aun(e))
          rhin(p1)  = rhin(p1)  - aneg
          rhout(p1) = rhout(p1) + apos
          rhin(p2)  = rhin(p2)  + apos
          rhout(p2) = rhout(p2) - aneg
        end do

        do inode=1,nb_nodes
          assert(vol(inode) .ne. 0._jprb )
          cp(inode) = ( Qmax(inode)-Q(inode) )*vol(inode)/( rhin(inode) * dt + eps )
          cn(inode) = ( Q(inode)-Qmin(inode) )*vol(inode)/( rhout(inode)* dt + eps )
        enddo

       ! limited antidiffusive  velocities:

        do LOOP_EDGES(e)
          p1 = edges(e,1)
          p2 = edges(e,2)
          if(aun(e) > 0.) then
            aun(e)=aun(e)*min(1._jprb,cp(p2),cn(p1))
          else
            aun(e)=aun(e)*min(1._jprb,cn(p2),cp(p1))
          endif
        enddo
      endif

      ! Compute fluxes from (limited) antidiffusive velocity
      do LOOP_EDGES(e)
        p1 = edges(e,1)
        p2 = edges(e,2)
        flux = aun(e)
        adv(p1) = adv(p1) + flux
        adv(p2) = adv(p2) - flux
      end do

      ! Update the unknowns in vertices
      do inode=1,nb_nodes
        Q(inode) = Q(inode) - adv(inode)/vol(inode) * dt
      end do

    end do ! other passes
    end associate
  end subroutine mpdata_gauge

  subroutine compute_gradient(Q,gradQ,Q_is_vector,geom)
    type(Geometry_class), intent(inout) :: geom
    real(kind=jprb), intent(in)    :: Q(:)
    real(kind=jprb), intent(inout) :: gradQ(:,:)
    logical, intent(in) :: Q_is_vector
    
    real(kind=jprb) :: Sx,Sy,avgQ
    integer :: iedge,p1,p2, e

    ! derivatives 
    gradQ(:,:) = 0

    do LOOP_EDGES(e)
      p1  = geom%edges(e,1)
      p2  = geom%edges(e,2)
      Sx  = geom%dual_normals(e,XX)
      Sy  = geom%dual_normals(e,YY)
      avgQ = ( Q(p1) + Q(p2) )*0.5_jprb
      gradQ(p1,XX) = gradQ(p1,XX) + Sx*avgQ
      gradQ(p2,XX) = gradQ(p2,XX) - Sx*avgQ
      gradQ(p1,YY) = gradQ(p1,YY) + Sy*avgQ
      gradQ(p2,YY) = gradQ(p2,YY) - Sy*avgQ
    end do

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    if (.not. Q_is_vector) then
      do iedge = 1,geom%nb_pole_edges
        e    = geom%pole_edges(iedge)
        p1   = geom%edges(e,1)
        p2   = geom%edges(e,2)
        Sy   = geom%dual_normals(e,YY)
        avgQ = ( Q(p1) + Q(p2) )*0.5_jprb
        ! correct for wrongly subtracting sy*avgQ in previous loop,
        ! instead of adding, because point at other side of pole
#ifdef  OLD_WAY
        gradQ(p1,YY) = gradQ(p1,YY) + Sy*avgQ
        gradQ(p2,YY) = gradQ(p2,YY) + Sy*avgQ
#else
        gradQ(p2,YY) = gradQ(p2,YY) + 2*Sy*avgQ 
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
  real(kind=jprb), parameter :: pi     = acos(-1.)


  real(kind=jprb) :: dt_forward = 20.
  integer :: iter = 0

contains
 
  subroutine set_time_step(dt)
    real(kind=jprb), intent(in) :: dt
    dt_forward = dt
  end subroutine set_time_step

  subroutine compute_metrics(geom)
    type(Geometry_class), intent(inout) :: geom
    real(kind=jprb) :: y, G, cos_y, sin_y
    integer :: inode
    
    write(0,*) "compute metrics"
    do inode=1,geom%nb_nodes
      y = geom%coordinates(inode,YY)
      cos_y = cos(y)
      sin_y = sin(y)
      geom%hx(inode) = radius*cos_y
      geom%hy(inode) = radius
      geom%dhxdy_over_G = -sin_y/(radius*max(eps,cos_y))
      G = geom%hx(inode) * geom%hy(inode)
      geom%dual_volumes(inode) = geom%dual_volumes(inode)*G
      geom%coriolis(inode) = f0 * sin_y
    enddo
    write(0,*) "compute metrics exit"
  end subroutine compute_metrics



  subroutine allocations(geom)
    type(Geometry_class), intent(inout) :: geom
    call create_state_fields(geom)
  end subroutine allocations



  subroutine propagate_state(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(Geometry_class), intent(inout), target :: geom
    real(kind=jprb) :: tmax, t0, dt_fwd
    tmax = geom%fields%time+dt
    do while (geom%fields%time < tmax)
      t0 = geom%fields%time
      dt_fwd = min( dt_forward, tmax-t0 )
      call step_forward(iter,dt_fwd,geom)

      write(0,'(A6,I8,A12,F9.1,A12,F8.1)') "step = ",iter, &
         & "  time = ",geom%fields%time, &
         & "  dt = ",dt_fwd
    end do
  end subroutine propagate_state




  subroutine step_forward(step,dt,geom)
    integer, intent(inout) :: step
    real(kind=jprb), intent(in) :: dt
    type(Geometry_class), intent(inout) :: geom

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
    type(Geometry_class), intent(inout) :: geom
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
    type(Geometry_class), intent(inout)      :: geom
    real(kind=jprb), dimension(:), pointer   :: D
    real(kind=jprb), dimension(:,:), pointer :: Q
    integer :: inode
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

    do inode=1,geom%nb_nodes
      x=geom%coordinates(inode,XX)
      y=geom%coordinates(inode,YY)
      if(x == 2._jprb*pi) x=0.
      cor=f0*sin(y)
      Q(inode,XX) =  radius*OM*cos(y)+radius*ZK*cos(IR*x) *(cos(y))**(IR-1)*(IR*(sin(y))**2-(cos(y))**2)
      Q(inode,YY) = -radius*ZK*IR*(cos(y))**(IR-1)*sin(y)*sin(IR*x)
      D(inode) = (ph0+radius**2*ATH(y)+radius**2*BTH(y)*cos(IR*x)+radius**2*CTH(y)*cos(2._jprb*IR*x)) &
          & /(grav)
      D(inode) = max(aaa0,D(inode))
      Q(inode,XX) = Q(inode,XX) * D(inode)
      Q(inode,YY) = Q(inode,YY) * D(inode)
      if(y == 0.5_jprb*pi) Q(inode,XX)=0.
      if(y ==-0.5_jprb*pi) Q(inode,XX)=0.
    end do
  end subroutine set_state_rossby_haurwitz



  subroutine backup_solution(geom)
    type(Geometry_class), intent(inout) :: geom
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
    type(Geometry_class), intent(inout) :: geom
    real(kind=jprb) :: Qx, Qy, Q0x, Q0y, Dmod, D0mod
    integer :: inode, iedge, e, p1, p2
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
    do inode=1,geom%nb_nodes
      Qx    = Q(inode,XX)
      Qy    = Q(inode,YY)
      Dmod  = max( eps, D(inode) )
      Q0x   = Q0(inode,XX)
      Q0y   = Q0(inode,YY)
      D0mod = max( eps, D0(inode) )
      
      ! this really computes V = G*contravariant_velocity, 
      ! with    G=hx*hy,
      !         physical_velocity = dotproduct( [hx,hy] , contravariant_velocity )
      ! V = (hx*hy) * [u/hx, v/hy] = [u*hy, v*hx]
      ! and hx = r*cos(y)  ,  hy = r
      ! and Q = [ D*u , D*v ]
      Vnodes(inode,XX) = ( 1.5_jprb*Qx/Dmod - 0.5_jprb*Q0x/D0mod )*geom%hy(inode)
      Vnodes(inode,YY) = ( 1.5_jprb*Qy/Dmod - 0.5_jprb*Q0y/D0mod )*geom%hx(inode)
    end do

    write(0,*) L2norm(Vnodes(:,XX)), L2norm(Vnodes(:,YY))

#ifndef OLD_WAY
    do e=1,geom%nb_edges
      p1 = geom%edges(e,1)
      p2 = geom%edges(e,2)
      Vedges(e,XX) = 0.5*(Vnodes(p1,XX)+Vnodes(p2,XX))
      Vedges(e,YY) = 0.5*(Vnodes(p1,YY)+Vnodes(p2,YY))
    enddo
    do iedge=1,geom%nb_pole_edges
      e = geom%pole_edges(iedge)
      Vedges(e,YY) = 0.
    enddo
#endif
  end subroutine compute_advective_velocities

  subroutine compute_forcing(geom)
    type(Geometry_class), intent(inout) :: geom
    integer :: inode
    real(kind=jprb) :: f, x, Qx, Qy, Dmod, dDdx, dDdy, vol, hx, hy
    real(kind=jprb), dimension(:),   pointer :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    real(kind=jprb) :: grad_D(geom%nb_nodes, 2)
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)

    call compute_gradient( D, grad_D, .False., geom )
    
    do inode=1,geom%nb_nodes
      Qx    = Q(inode,XX)
      Qy    = Q(inode,YY)
      Dmod  = max( eps, D(inode) )
      vol   = geom%dual_volumes(inode)
      hx    = geom%hx(inode)
      hy    = geom%hy(inode)
      f     = geom%coriolis(inode)
      R(inode,XX) = -grav*D(inode)*grad_D(inode,XX)*hy/vol + f*Qy - geom%dhxdy_over_G(inode)*Qx*Qy/Dmod
      R(inode,YY) = -grav*D(inode)*grad_D(inode,YY)*hx/vol - f*Qx + geom%dhxdy_over_G(inode)*Qx*Qx/Dmod
    end do

  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(Geometry_class), intent(inout) :: geom
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    integer :: inode, nb_nodes
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)
    do inode=1,geom%nb_nodes
      Q(inode,XX) = Q(inode,XX) + 0.5_jprb*dt*R(inode,XX)
      Q(inode,YY) = Q(inode,YY) + 0.5_jprb*dt*R(inode,YY)
    end do
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(Geometry_class), intent(inout) :: geom
    integer :: inode, m
    real(kind=jprb) :: f, Qx, Qy, Rx, Ry, Dmod, dDdx, dDdy, hx, hy
    real(kind=jprb) :: Qx_adv, Qy_adv, Rx_exp, Ry_exp, vol

    real(kind=jprb), dimension(:),   pointer :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    real(kind=jprb) :: grad_D(geom%nb_nodes, 2)

    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    call compute_gradient( D, grad_D, .False., geom )
    
    do inode=1,geom%nb_nodes
      Qx    = Q(inode,XX)
      Qy    = Q(inode,YY)
      Dmod  = max( eps, D(inode) )
      vol   = geom%dual_volumes(inode)
      hx    = geom%hx(inode)
      hy    = geom%hy(inode)
      f     = geom%coriolis(inode)

      Rx_exp = -grav*D(inode)*grad_D(inode,XX)*hy/vol
      Ry_exp = -grav*D(inode)*grad_D(inode,YY)*hx/vol

      Qx_adv = Qx
      Qy_adv = Qy

      do m=1,3 ! Three iterations at most is enough to converge
        Rx = Rx_exp + f*Qy - geom%dhxdy_over_G(inode)*Qx*Qy/Dmod
        Ry = Ry_exp - f*Qx + geom%dhxdy_over_G(inode)*Qx*Qx/Dmod
        Qx = Qx_adv + 0.5_jprb*dt*Rx
        Qy = Qy_adv + 0.5_jprb*dt*Ry
      end do
      Q(inode,XX) = Qx
      Q(inode,YY) = Qy
      R(inode,XX) = Rx_exp + f*Qy - geom%dhxdy_over_G(inode)*Qx*Qy/Dmod
      R(inode,YY) = Ry_exp - f*Qx + geom%dhxdy_over_G(inode)*Qx*Qx/Dmod
    end do
  end subroutine implicit_solve



  subroutine advect_solution(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(Geometry_class), intent(inout) :: geom
    real(kind=jprb), dimension(:),   pointer :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, V
    integer :: inode, p1, p2
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    V => vector_field("advective_velocity",geom)

    !    mpdata_gauge( time, variable, velocity, order, limit,  is_vector, geom )
    call mpdata_gauge( dt,   D,        V,        2,     .True., .False.,   geom )
    call mpdata_gauge( dt,   Q(:,XX),  V,        2,     .True., .True. ,   geom )
    call mpdata_gauge( dt,   Q(:,YY),  V,        2,     .True., .True. ,   geom )

    ! remove noise from periodic boundary.. why is it even there in the mesh?
    do inode=1,geom%nb_ghost_nodes
      p1 = geom%ghost_nodes(inode,1)
      p2 = geom%ghost_nodes(inode,2)
      D(p2)    = D(p1)
      Q(p2,XX) = Q(p1,XX)
      Q(p2,YY) = Q(p1,YY)
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
  use datastruct_module,  only: Geometry_class
  use shallow_water_module, only: &
    & compute_metrics, &
    & allocations, &
    & set_state_rossby_haurwitz, &
    & set_time_step, &
    & propagate_state
  implicit none

  ! Configuration parameters
  real(kind=jprb) :: dt = 20.             ! solver time-step
  integer         :: nb_steps = 6        ! Number of propagations
  integer         :: hours_per_step = 1   !Propagation time

  ! Declarations
  type(Geometry_class) :: g
  real(kind=jprb), parameter :: hours = 3600.     ! in seconds
  real(kind=jprb), parameter :: days  = 24.*hours ! in seconds
  integer :: step = 0
  character(len=1024) :: filename
  integer :: clck_counts_start, clck_counts_beg, clck_counts_end, clck_rate
 
  ! Execution
  call read_joanna("data/meshvol.d", g)
  call compute_metrics(g)
  call allocations(g)
  call set_state_rossby_haurwitz(g)
  call set_time_step( dt )

  call write_gmsh_mesh(g%internal_mesh,"data/mesh.msh")

  write (filename, "(A,I2.2,A)") "data/fields",step,".msh"
  call write_gmsh_state(g%fields,filename)

  write (filename, "(A,I2.2,A)") "data/depth",step,".grib"
  call write_grib(g,filename)

  call system_clock ( clck_counts_start, clck_rate )
  do step=1,nb_steps 

    call system_clock ( clck_counts_beg, clck_rate )
  
    call propagate_state( hours_per_step*hours, g)
    !call propagate_state( 20._jprb, g)
  
    call system_clock ( clck_counts_end, clck_rate )
    write (0, *) "propagated to ",step*hours_per_step," hours.   time = ",(clck_counts_end - clck_counts_beg) / real(clck_rate),&
     &  "walltime = ",(clck_counts_end - clck_counts_start) / real(clck_rate), new_line('A')
    
    write (filename, "(A,I2.2,A)") "data/fields",step,".msh"
    call write_gmsh_state(g%fields,filename)

    write (filename, "(A,I2.2,A)") "data/depth",step,".grib"
    call write_grib(g,filename)
  end do
  
end program shallow_water
