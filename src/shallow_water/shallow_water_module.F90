! ===================================================================
! shallow_water_module
! --------------------
! This module contains subroutines to solve the shallow water eqs
! on a sphere
! - set_time_step : To set the solver's time step
! - setup_shallow_water : To modify dual-volumes to map to the sphere
!                         and create necessary fields
! - set_state_rossby_haurwits : To initialise the state with
!                               Rossby Haurwitz waves
! - propagate_state : To propagate the state with a given time step
! ===================================================================
module shallow_water_module
  use parallel_module
  use common_module
  use datastruct_module
  use mpdata_module

  implicit none
  private
  public :: setup_shallow_water 
  public :: propagate_state
  public :: set_state_rossby_haurwitz
  public :: set_state_zonal_flow
  public :: set_topography_mountain
  public :: set_time_step

  integer, parameter, public :: EQS_MOMENTUM = 1
  integer, parameter, public :: EQS_VELOCITY = 2

  real(kind=jprw), parameter :: eps    = 1.e-6
  real(kind=jprw), parameter :: radius = 6371.22e+03
!  real(kind=jprw), parameter :: radius = 63.7122e+03

  real(kind=jprw), parameter :: f0     = 1.4584e-04 !coriolis parameter (=2xearth's omega)
  real(kind=jprw), parameter :: grav   = 9.80616
  real(kind=jprw), parameter :: pi     = acos(-1._jprw)
  real(kind=jprw), parameter :: D_tres = 1._jprw

  real(kind=jprw) :: dt_forward = 20.
  integer :: iter = 0
  integer :: eqs_type
contains
 
  subroutine set_time_step(dt)
    real(kind=jprw), intent(in) :: dt
    dt_forward = dt
  end subroutine set_time_step

  subroutine setup_shallow_water(eqs_type_,dstruct)
    integer, intent(in) :: eqs_type_
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw) :: y, cos_y, sin_y
    real(kind=jprw), pointer :: coords(:,:), vol(:), hx(:), hy(:), dhxdy_over_G(:), pole_bc(:), G(:)
    integer :: jnode, jedge, iedge

    eqs_type = eqs_type_

    call create_field_in_nodes_2d("jacobian",1,dstruct)
    call create_field_in_nodes_2d("hx",1,dstruct)
    call create_field_in_nodes_2d("hy",1,dstruct)
    call create_field_in_nodes_2d("dhxdy_over_G",1,dstruct)
    call create_field_in_edges_2d("pole_bc",1,dstruct)


    coords => vector_field_2d("coordinates",dstruct)
    vol    => scalar_field_2d("dual_volumes",dstruct)
    hx    => scalar_field_2d("hx",dstruct)
    hy    => scalar_field_2d("hy",dstruct)
    G     => scalar_field_2d("jacobian",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)

    !dir$ ivdep
    do jnode=1,dstruct%nb_nodes
      y = coords(YY,jnode)
      cos_y = cos(y)
      sin_y = sin(y)
      hx(jnode) = radius*cos_y
      hy(jnode) = radius
      G(jnode) = hx(jnode)*hy(jnode)
      vol(jnode) = vol(jnode)*G(jnode)
      dhxdy_over_G(jnode) = - sin_y/(radius*max(eps,cos_y))
    enddo

    pole_bc(:) = 1.
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      pole_bc(iedge) = -1.
    end do

    call create_field_in_nodes_2d("coriolis",1,dstruct)
    call create_field_in_nodes_2d("depth",1,dstruct)
    call create_field_in_nodes_2d("velocity",2,dstruct)
    call create_field_in_nodes_2d("momentum",2,dstruct)
    call create_field_in_nodes_2d("forcing",2,dstruct)
    call create_field_in_nodes_2d("depth_backup",1,dstruct)
    call create_field_in_nodes_2d("velocity_backup",2,dstruct)
    call create_field_in_nodes_2d("momentum_backup",2,dstruct)
    call create_field_in_nodes_2d("initial_velocity",2,dstruct)
    call create_field_in_edges_2d("advective_velocity",2,dstruct)
    call create_field_in_nodes_2d("depth_ratio",1,dstruct) ! old/new
    call create_field_in_nodes_2d("divV",1,dstruct)
    
    call create_field_in_nodes_2d("topography",1,dstruct)
    call create_field_in_nodes_2d("height",1,dstruct)

  end subroutine setup_shallow_water

  subroutine propagate_state(dt,order,scheme,dstruct)
    real(kind=jprw), intent(in) :: dt
    integer, intent(in) :: order
    integer, intent(in) :: scheme
    type(DataStructure_type), intent(inout), target :: dstruct
    real(kind=jprw) :: tend, t0, dt_fwd, tstart
    character(len=200) :: step_info
    real(kind=jprw), pointer :: D(:)

    D => scalar_field_2d("depth",dstruct)
    tstart   = dstruct%time
    tend     = dstruct%time+dt
    
    call log_info(" ")
    write(log_str,'(A,I10,A,I10)') "Propagating from time ",int(tstart)," to time ",int(tend); call log_info()
    do while (dstruct%time < tend)
      t0 = dstruct%time
      dt_fwd = min( dt_forward, tend-t0 )
      call step_forward(iter,dt_fwd,order,scheme,dstruct)

      if( log_level <= LOG_LEVEL_INFO ) then
        if( myproc .eq. 0 ) then
          call progress_bar(dstruct%time,tstart,tend)
        end if
      else
        write(step_info,'(A6,I8,A12,F9.1,A12,F8.1,A12,E20.13)') "step = ",iter, &
          & "  time = ",dstruct%time, &
          & "  dt = ",dt_fwd, "Norm ",L2norm(D)
        CALL LOG_INFO( STEP_INFO )
      end if
    end do
    dstruct%time_step = dstruct%time_step + 1

  end subroutine propagate_state




  subroutine step_forward(step,dt,order,scheme,dstruct)
    integer, intent(inout) :: step
    real(kind=jprw), intent(in) :: dt
    integer, intent(in) :: order
    integer, intent(in) :: scheme
    type(DataStructure_type), intent(inout) :: dstruct

    call backup_solution(dstruct)

    if (step == 0) then ! Pre-compute forcing

      call compute_forcing(dstruct)
      call compute_advective_velocities(dt,dstruct,"linear_advection")

    end if

    call add_forcing_to_solution(dt,dstruct)

    call advect_solution(dt,order,scheme,dstruct)

    call implicit_solve(dt,dstruct)

    !call filter_solution(dstruct)

    call compute_advective_velocities(dt,dstruct,"linear_advection")

    dstruct%time = dstruct%time+dt
    step = step+1

  end subroutine step_forward



  subroutine set_state_rossby_haurwitz(dstruct)
    type(DataStructure_type), intent(inout)      :: dstruct
    real(kind=jprw), dimension(:), pointer   :: D, cor, H0, H
    real(kind=jprw), dimension(:,:), pointer :: U, Q, coords
    integer :: jnode, ir
    real(kind=jprw) :: zk,om,ph0,x,y, sin_y, cos_y

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3

    coords => vector_field_2d("coordinates",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    Q => vector_field_2d("momentum",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_2d("height",dstruct)


    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y,sin_y,cos_y)
    do jnode=1,dstruct%nb_nodes
      x=coords(XX,jnode)
      y=coords(YY,jnode)
      sin_y = sin(y)
      cos_y = cos(y)
      cor(jnode) = f0*sin_y
      if(x == 2._jprw*pi) x=0.
      U(XX,jnode) =  radius*om*cos_y+radius*zk*cos(ir*x) * cos_y**(ir-1._jprw) &
        &            * (ir*(sin_y)**2-cos_y**2)
      U(YY,jnode) = -radius*zk*ir*cos_y**(ir-1._jprw)*sin_y*sin(ir*x)
      H(jnode) = (ph0+radius**2*fa(y)+radius**2*fb(y)*cos(ir*x) &
        &        +radius**2*fc(y)*cos(2._jprw*ir*x)) / grav
      D(jnode) = max(0._jprw,H(jnode) - H0(jnode))
      if(y == 0.5_jprw*pi) U(XX,jnode)=0.
      if(y ==-0.5_jprw*pi) U(XX,jnode)=0.
      H(jnode) = H0(jnode) + D(jnode)
      Q(:,jnode) = D(jnode) * U(:,jnode)
    end do
    !$OMP END PARALLEL DO

    contains 
    ! Helper functions

      real(kind=jprw) function fa(th)
        real(kind=jprw), intent(in) :: th
        fa = om*0.5*(f0+om)*(cos(th))**2 &
          & +0.25*zk**2*(cos(th))**(2*ir)*( (ir+1)*(cos(th))**2 &
          & +(2._jprw*ir**2-ir-2._jprw)-2._jprw*ir**2/(cos(th))**2 )
      end function fa

      real(kind=jprw) function fb(th)
        real(kind=jprw), intent(in) :: th
        fb = (f0+2._jprw*om)*zk/((ir+1._jprw)*(ir+2._jprw))*(cos(th))**ir &
          & *( (ir**2+2._jprw*ir+2._jprw)-((ir+1._jprw)*cos(th))**2 )
      end function fb

      real(kind=jprw) function fc(th)
        real(kind=jprw), intent(in) :: th
        fc = 0.25*zk**2*(cos(th))**(2*ir)*((ir+1._jprw)*(cos(th))**2 -(ir+2._jprw))  
      end function fc

  end subroutine set_state_rossby_haurwitz


  subroutine set_state_zonal_flow(dstruct)
    type(DataStructure_type), intent(inout)  :: dstruct
    real(kind=jprw), dimension(:), pointer   :: D, cor, H0, H
    real(kind=jprw), dimension(:,:), pointer :: U, Uinit, Q, coords
    integer :: jnode
    real(kind=jprw) :: x,y
    real(kind=jprw), parameter :: USCAL = 10.
    real(kind=jprw), parameter :: H00 = grav * 8e3
    real(kind=jprw), parameter :: pvel = USCAL/radius
    real(kind=jprw), parameter :: beta = pi/4._jprw


    coords => vector_field_2d("coordinates",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    Q => vector_field_2d("momentum",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_2d("height",dstruct)
    Uinit => vector_field_2d("initial_velocity",dstruct)

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
    do jnode=1,dstruct%nb_nodes
      x=coords(XX,jnode)
      y=coords(YY,jnode)
      cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
      H(jnode)     = (H00-radius**2*(f0+pvel)*0.5*pvel &
                   & *(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2)/grav
      D(jnode)     = max(0., H(jnode) - H0(jnode) )
      U(XX,jnode)  =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y)
      U(YY,jnode)  = -pvel*sin(x)*sin(beta)*radius
      if ( D(jnode) < D_tres ) then
        D(jnode) = 0.
        U(:,jnode) = 0.
      end if
      H(jnode) = H0(jnode) + D(jnode)
      Q(:,jnode) = D(jnode) * U(:,jnode)
      Uinit(:,jnode) = U(:,jnode)
    end do
    !$OMP END PARALLEL DO
    call spectral_filter(Uinit,1._jprb,dstruct)
    call compute_forcing(dstruct)


  end subroutine set_state_zonal_flow


  subroutine set_topography_mountain(amplitude,dstruct)
    real(kind=jprw), intent(in) :: amplitude ! amplitude of hill
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:), pointer :: H0
    real(kind=jprw), dimension(:,:), pointer :: coords
    real(kind=jprw) :: rad = 2.*pi/18. ! radius of hill
    real(kind=jprw) :: xcent = 3.*pi/2.  ! centre of hill
    real(kind=jprw) :: ycent = pi/6.*1.
    real(kind=jprw) :: gamm = 1. ! slope of hill
    real(kind=jprw) :: dist, xlon, ylat
    integer :: jnode

    H0 => scalar_field_2d("topography",dstruct)
    coords => vector_field_2d("coordinates",dstruct)
    
    do jnode=1,dstruct%nb_nodes
      xlon = coords(XX,jnode)
      ylat = coords(YY,jnode)

      dist = 2.*sqrt( (cos(ylat)*sin( (xlon-xcent)/2 ) )**2 &
        &     + sin((ylat-ycent)/2)**2 )

      ! This is a cosine hill
      H0(jnode) = 0.
      if (dist < rad) H0(jnode) = 0.5*(1.+cos(pi*dist/rad))*amplitude

      ! This is a cone
      !H0(jnode) = max(0.,amplitude * (1.-gamm*dist/rad))

      ! This is a cylinder
      !H0(jnode) = 0.
      !if (dist < rad) H0(jnode) = amplitude

    end do
  end subroutine set_topography_mountain


  subroutine backup_solution(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:),   pointer :: D, D0
    real(kind=jprw), dimension(:,:), pointer :: U, U0
    real(kind=jprw), dimension(:,:), pointer :: Q, Q0
    integer :: jnode
    
    D  => scalar_field_2d("depth",dstruct)
    D0 => scalar_field_2d("depth_backup",dstruct)
    U  => vector_field_2d("velocity",dstruct)
    U0 => vector_field_2d("velocity_backup",dstruct)
    Q  => vector_field_2d("momentum",dstruct)
    Q0 => vector_field_2d("momentum_backup",dstruct)

    select case (eqs_type)

      case (EQS_MOMENTUM)
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
        do jnode=1,dstruct%nb_nodes
          D0(jnode)   = D(jnode)
          Q0(:,jnode) = Q(:,jnode)
        end do
        !$OMP END PARALLEL DO

      case (EQS_VELOCITY)

        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
        do jnode=1,dstruct%nb_nodes
          D0(jnode)   = D(jnode)
          U0(:,jnode) = U(:,jnode)
        end do
       !$OMP END PARALLEL DO

    end select

  end subroutine backup_solution



  subroutine compute_advective_velocities(dt,dstruct,option)
    ! this really computes V = G*contravariant_velocity, 
    ! with    G=hx*hy,
    !         physical_velocity = dotproduct( [hx,hy] , contravariant_velocity )
    ! V = (hx*hy) * [u/hx, v/hy] = [u*hy, v*hx]
    ! and hx = r*cos(y)  ,  hy = r
    ! and Q = [ D*u , D*v ]
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    character(len=*), intent(in), optional :: option
    real(kind=jprw) :: Ux, Uy, U0x, U0y, Vx, Vy, dVxdx, dVxdy, dVydx, dVydy
    real(kind=jprw) :: Dpos, D0pos
    integer :: jnode, jedge, iedge, ip1, ip2
    real(kind=jprw), dimension(:),   pointer :: D, D0, hx, hy, vol
    real(kind=jprw), dimension(:,:), pointer :: U, U0, R, Vedges, coords, Q, Q0, S, Uinit
    real(kind=jprw) :: Vnodes(2,dstruct%nb_nodes), grad_Vnodes(4,dstruct%nb_nodes)
    real(kind=jprw), parameter :: LAGR = 0.

    coords => vector_field_2d("coordinates",dstruct)
    Vedges => vector_field_2d("advective_velocity",dstruct)
    D      => scalar_field_2d("depth",dstruct)
    D0     => scalar_field_2d("depth_backup",dstruct)
    Q      => vector_field_2d("momentum",dstruct)
    Q0     => vector_field_2d("momentum_backup",dstruct)
    U      => vector_field_2d("velocity",dstruct)
    Uinit  => vector_field_2d("initial_velocity",dstruct)
    U0     => vector_field_2d("velocity_backup",dstruct)
    hx     => scalar_field_2d("hx",dstruct)
    hy     => scalar_field_2d("hy",dstruct)
    R      => vector_field_2d("forcing",dstruct)
    vol    => scalar_field_2d("dual_volumes",dstruct)
    S      => vector_field_2d("dual_normals",dstruct)

    if( option .eq. "advect") then
      select case (eqs_type)

        case (EQS_MOMENTUM)
          !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Dpos)
          do jnode=1,dstruct%nb_nodes
            Dpos = max(eps, D(jnode))
            Vnodes(:,jnode)=Q(:,jnode)/Dpos + LAGR*R(:,jnode)/Dpos
            U(:,jnode) = Q(:,jnode)/Dpos
          end do
          !$OMP END PARALLEL DO
        case (EQS_VELOCITY)
          !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
          do jnode=1,dstruct%nb_nodes
            Vnodes(:,jnode)=(U(:,jnode)+0.5*dt*R(:,jnode)) 
          end do
          !$OMP END PARALLEL DO
      end select

      call compute_gradient_tensor( Vnodes, grad_Vnodes, dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Vx,Vy,dVxdx,dVxdy,dVydx,dVydy)
      do jnode=1,dstruct%nb_nodes
        Ux = U(XX,jnode)
        Uy = U(YY,jnode)
        Vx = Ux!/hx(jnode)  !bad results when uncommented
        Vy = Uy!/hy(jnode)
        dVxdx = grad_Vnodes(XXDXX,jnode)*hy(jnode)/vol(jnode)
        dVxdy = grad_Vnodes(XXDYY,jnode)*hx(jnode)/vol(jnode)
        dVydx = grad_Vnodes(YYDXX,jnode)*hy(jnode)/vol(jnode)    
        dVydy = grad_Vnodes(YYDYY,jnode)*hx(jnode)/vol(jnode)
        Dpos = max(eps,D(jnode))
        Vnodes(XX,jnode) = ( Ux - 0.5*dt*(Vx*dVxdx+Vy*dVxdy - (1.-LAGR)*R(XX,jnode)/Dpos ) ) * hy(jnode)
        Vnodes(YY,jnode) = ( Uy - 0.5*dt*(Vx*dVydx+Vy*dVydy - (1.-LAGR)*R(YY,jnode)/Dpos ) ) * hx(jnode)
      end do
      !$OMP END PARALLEL DO

      call halo_exchange( Vnodes, dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge=1,dstruct%nb_edges
        ip1 = dstruct%edges(1,jedge)
        ip2 = dstruct%edges(2,jedge)
        if (DENSITY_WEIGHTING) then
          Vedges(:,jedge) = (D(ip1)*Vnodes(:,ip1)+D(ip2)*Vnodes(:,ip2))/(D(ip1)+D(ip2)+eps)
        else
          Vedges(:,jedge) = (Vnodes(:,ip1)+Vnodes(:,ip2))*0.5
        end if
      end do
      !$OMP END PARALLEL DO

    else if( option .eq. "extrapolate") then
      select case (eqs_type)

        case (EQS_MOMENTUM)
          !dir$ ivdep
          !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Dpos,D0pos,Ux,Uy,U0x,U0y)
          do jnode=1,dstruct%nb_nodes
            Dpos  = max(eps, D(jnode))
            D0pos  = max(eps, D0(jnode))
            Ux    = Q(XX,jnode)/Dpos
            Uy    = Q(YY,jnode)/Dpos
            U0x   = Q0(XX,jnode)/D0pos
            U0y   = Q0(YY,jnode)/D0pos
            Vnodes(XX,jnode) = ( 1.5_jprw*Ux - 0.5_jprw*U0x ) * hy(jnode)
            Vnodes(YY,jnode) = ( 1.5_jprw*Uy - 0.5_jprw*U0y ) * hx(jnode)
          end do
          !$OMP END PARALLEL DO

        case (EQS_VELOCITY)
          !dir$ ivdep
          !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy,U0x,U0y)
          do jnode=1,dstruct%nb_nodes
            Ux    = U(XX,jnode)
            Uy    = U(YY,jnode)
            U0x   = U0(XX,jnode)
            U0y   = U0(YY,jnode)
            Vnodes(XX,jnode) = ( 1.5_jprw*Ux - 0.5_jprw*U0x ) * hy(jnode)
            Vnodes(YY,jnode) = ( 1.5_jprw*Uy - 0.5_jprw*U0y ) * hx(jnode)
          end do
          !$OMP END PARALLEL DO
      end select

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge=1,dstruct%nb_edges
        ip1 = dstruct%edges(1,jedge)
        ip2 = dstruct%edges(2,jedge)
        if (DENSITY_WEIGHTING) then
          Vedges(:,jedge) = (D(ip1)*Vnodes(:,ip1)+D(ip2)*Vnodes(:,ip2))/(D(ip1)+D(ip2)+eps)
        else
          Vedges(:,jedge) = (Vnodes(:,ip1)+Vnodes(:,ip2))*0.5
        end if
      enddo
      !$OMP END PARALLEL DO

    else if( option .eq. "linear_advection") then
      do jnode=1,dstruct%nb_nodes
        Vnodes(XX,jnode) = Uinit(XX,jnode) * hy(jnode)
        Vnodes(YY,jnode) = Uinit(YY,jnode) * hx(jnode)
      end do
      do jedge=1,dstruct%nb_edges
        ip1 = dstruct%edges(1,jedge)
        ip2 = dstruct%edges(2,jedge)
        if (DENSITY_WEIGHTING) then
          Vedges(:,jedge) = (D(ip1)*Vnodes(:,ip1)+D(ip2)*Vnodes(:,ip2))/(D(ip1)+D(ip2)+eps)
        else
          Vedges(:,jedge) = (Vnodes(:,ip1)+Vnodes(:,ip2))*0.5
        end if
      enddo

    end if

   
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      Vedges(YY,iedge) = 0.
    enddo

  end subroutine compute_advective_velocities

  subroutine compute_forcing(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode
    real(kind=jprw) :: Ux, Uy, Qx, Qy, Dpos
    real(kind=jprw), dimension(:),   pointer :: H, H0, D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:), pointer :: U, Q, R
    real(kind=jprw) :: grad_H(2,dstruct%nb_nodes)
    vol => scalar_field_2d("dual_volumes",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    H => scalar_field_2d("height",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    Q => vector_field_2d("momentum",dstruct)
    R => vector_field_2d("forcing",dstruct)

    H(:) = H0(:) + D(:)
    call compute_gradient( H, grad_H, dstruct )
    call halo_exchange(grad_H,dstruct)

    select case (eqs_type)
      case (EQS_MOMENTUM)
        !dir$ ivdep
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
        do jnode=1,dstruct%nb_nodes
          Qx    = Q(XX,jnode)
          Qy    = Q(YY,jnode)
          Dpos  = max(eps, D(jnode))
          R(XX,jnode) = -grav*D(jnode)*grad_H(XX,jnode)*hy(jnode)/vol(jnode) &
            &           + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dpos
          R(YY,jnode) = -grav*D(jnode)*grad_H(YY,jnode)*hx(jnode)/vol(jnode) &
            &           - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dpos
        end do
        !$OMP END PARALLEL DO

      case (EQS_VELOCITY)
        !dir$ ivdep
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
        do jnode=1,dstruct%nb_nodes
          Ux    = U(XX,jnode)
          Uy    = U(YY,jnode)
          R(XX,jnode) = -grav*grad_H(XX,jnode)*hy(jnode)/vol(jnode) &
            &           + cor(jnode)*Uy - dhxdy_over_G(jnode)*Ux*Uy
          R(YY,jnode) = -grav*grad_H(YY,jnode)*hx(jnode)/vol(jnode) &
            &           - cor(jnode)*Ux + dhxdy_over_G(jnode)*Ux*Ux
        end do
        !$OMP END PARALLEL DO
    end select

  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: solution(:,:) , forcing(:,:)
    integer :: jnode

    forcing  => vector_field_2d("forcing",dstruct)

    select case (eqs_type)
      case (EQS_MOMENTUM)
        solution => vector_field_2d("momentum",dstruct)
      case (EQS_VELOCITY)
        solution => vector_field_2d("velocity",dstruct)
    end select

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      solution(:,jnode) = ( solution(:,jnode) + 0.5_jprw*dt*forcing(:,jnode) )
    end do
    !$OMP END PARALLEL DO
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, m
    real(kind=jprw) :: Ux, Uy, Qx, Qy, Rx, Ry, Dpos
    real(kind=jprw) :: Ux_adv, Uy_adv, Qx_adv, Qy_adv, Rx_exp, Ry_exp

    real(kind=jprw), dimension(:),   pointer :: H, H0, D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:), pointer :: U, Q, R
    real(kind=jprw) :: grad_H(2,dstruct%nb_nodes)

    vol => scalar_field_2d("dual_volumes",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_2d("height",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    Q => vector_field_2d("momentum",dstruct)
    R => vector_field_2d("forcing",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    H(:) = H0(:) + D(:)
    call compute_gradient( H, grad_H, dstruct )
    call halo_exchange( grad_H, dstruct )

    select case (eqs_type)

      case (EQS_MOMENTUM)
        !dir$ ivdep
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Qx,Qy,Dpos,Rx_exp,Ry_exp,Qx_adv,Qy_adv,m,Rx,Ry)
        do jnode=1,dstruct%nb_nodes
          Dpos  = D(jnode)

          Qx    = Q(XX,jnode)
          Qy    = Q(YY,jnode)

          Qx_adv = Qx
          Qy_adv = Qy

          Rx_exp = -grav*D(jnode)*grad_H(XX,jnode)*hy(jnode)/vol(jnode)
          Ry_exp = -grav*D(jnode)*grad_H(YY,jnode)*hx(jnode)/vol(jnode)

          if (D(jnode) > D_tres) then
            do m=1,3 ! Three iterations at most is enough to converge
              Rx = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/D(jnode)
              Ry = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/D(jnode)
              Qx = Qx_adv + 0.5_jprw*dt*Rx
              Qy = Qy_adv + 0.5_jprw*dt*Ry
            end do
            Q(XX,jnode) = Qx
            Q(YY,jnode) = Qy
            R(XX,jnode) = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/D(jnode)
            R(YY,jnode) = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/D(jnode)
            U(:,jnode)  = Q(:,jnode) / D(jnode)
          else
            R(:,jnode) = 0.
            U(:,jnode) = 0.
            Q(:,jnode) = 0.
          end if

          if ( jnode == probe ) then
            write(log_str,*) "Q_next ", jnode, Q(:,jnode); call log_debug()
            write(log_str,*) "U_next ", jnode, U(:,jnode); call log_debug()
            write(log_str,*) "R_exp  ", jnode, Rx_exp, Ry_exp; call log_debug()
            write(log_str,*) "R      ", jnode, R(:,jnode); call log_debug()
            write(log_str,*) "gradH  ", jnode, grad_H(:,jnode); call log_debug()
          end if

        end do
        !$OMP END PARALLEL DO



      case (EQS_VELOCITY)
        !dir$ ivdep
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy,Rx_exp,Ry_exp,Ux_adv,Uy_adv,m,Rx,Ry)
        do jnode=1,dstruct%nb_nodes
          Ux    = U(XX,jnode)
          Uy    = U(YY,jnode)

          Rx_exp = -grav*grad_H(XX,jnode)*hy(jnode)/vol(jnode)
          Ry_exp = -grav*grad_H(YY,jnode)*hx(jnode)/vol(jnode)

          Ux_adv = Ux
          Uy_adv = Uy

          do m=1,3 ! Three iterations at most is enough to converge
            Rx = Rx_exp + cor(jnode)*Uy - dhxdy_over_G(jnode)*Ux*Uy
            Ry = Ry_exp - cor(jnode)*Ux + dhxdy_over_G(jnode)*Ux*Ux
            Ux = Ux_adv + 0.5_jprw*dt*Rx
            Uy = Uy_adv + 0.5_jprw*dt*Ry
          end do
          U(XX,jnode) = Ux
          U(YY,jnode) = Uy
          R(XX,jnode) = Rx_exp + cor(jnode)*Uy - dhxdy_over_G(jnode)*Ux*Uy
          R(YY,jnode) = Ry_exp - cor(jnode)*Ux + dhxdy_over_G(jnode)*Ux*Ux
          Q(:,jnode) = U(:,jnode) * D(jnode)

          if (D(jnode) < D_tres) then
            U(:,jnode) = 0.
            R(:,jnode) = 0.
            Q(:,jnode) = 0.
          end if
        end do
        !$OMP END PARALLEL DO
    end select


  end subroutine implicit_solve

  subroutine filter_solution(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: D(:), U(:,:), Q(:,:)
    integer :: jnode
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    Q => vector_field_2d("momentum",dstruct)

    !call spectral_filter(D,0.67_jprw,dstruct)
    !call spectral_filter(Q,0.67_jprw,dstruct)

  end subroutine filter_solution


  subroutine advect_solution(dt,order,scheme,dstruct)
    real(kind=jprw), intent(in) :: dt
    integer, intent(in) :: order
    integer, intent(in) :: scheme
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:),   pointer :: D, D0, DR
    real(kind=jprw), dimension(:,:), pointer :: U, Q, V
    real(kind=jprw) :: VDS(dstruct%nb_edges)
    integer :: jnode

    real(kind=jprw), parameter :: limit = 1.
    
    D => scalar_field_2d("depth",dstruct)
    D0 => scalar_field_2d("depth_backup",dstruct)
    U => vector_field_2d("velocity",dstruct)
    Q => vector_field_2d("momentum",dstruct)
    V => vector_field_2d("advective_velocity",dstruct)
    DR => scalar_field_2d("depth_ratio",dstruct)
   
    !    mpdata_D( scheme, time, variable, velocity, VDS,  order, limit,   dstruct )
    call mpdata_D( scheme, dt,   D,        V,        VDS,  order, limit, dstruct )

    select case (eqs_type)

      case (EQS_MOMENTUM)
        !    mpdata_Q( scheme, time, variable, V,  order, limit, dstruct )
        call mpdata_Q( scheme, dt,   Q,        V,  order, limit, dstruct )

      case (EQS_VELOCITY)
        ! compute ratio
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
        do jnode=1,dstruct%nb_nodes
          DR(jnode) = D0(jnode) / max( D(jnode), eps )
        end do
        !$OMP END PARALLEL DO
        !    mpdata_gauge_U( time, variable, VDS, DR, D0,  order, limit,  dstruct )
        call mpdata_gauge_U( dt,   U,        VDS, DR, D0,  order, limit, dstruct )

    end select

  end subroutine advect_solution

end module shallow_water_module
