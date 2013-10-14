! =====================================================================
! mpdata_module
! This module contains strictly algorithmic subroutines
! - mpdata_gauge       : Non-oscillatory variable sign MPDATA advection
! - compute_gradient   : Compute gradient of a scalar array
! =====================================================================
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
    real(kind=jprb), pointer :: vol(:), S(:,:)

    associate( nb_nodes          => geom%nb_nodes, &
      &        nb_edges          => geom%nb_edges, &
      &        edges             => geom%edges )
  
    vol => scalar_field("dual_volumes",geom)
    S   => vector_field("dual_normals",geom)

    ! non-oscillatory option
    if( limited .and. (order .ge. 2) ) then
      call compute_Qmax_and_Qmin()
    end if

    ! 1. First pass
    ! -------------
    
    ! Compute the normal velocity in faces, and advection in vertices
    adv(:) = 0.

    do jedge = 1,nb_edges
      ip1 = edges(jedge,1)
      ip2 = edges(jedge,2)
      Sx = S(jedge,XX)
      Sy = S(jedge,YY)
      Vx = V(jedge,XX)
      Vy = V(jedge,YY)
      aun(jedge) = Vx*Sx + Vy*Sy
      apos = max(0._jprb,aun(jedge))
      aneg = min(0._jprb,aun(jedge))
      flux = Q(ip1)*apos + Q(ip2)*aneg
      adv(ip1) = adv(ip1) + flux
      adv(ip2) = adv(ip2) - flux
    end do

    ! Update the unknowns in vertices
    do jnode=1,nb_nodes
      Q(jnode) = Q(jnode) - adv(jnode)/vol(jnode) * dt
    end do


    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order
      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector, geom)

      ! Compute antidiffusive normal velocity in faces

      do jedge = 1,nb_edges
        ip1 = edges(jedge,1)
        ip2 = edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(ip1) + vol(ip2)
        dQdx = (gradQ(ip1,XX)+gradQ(ip2,XX)) / volume_of_two_cells
        dQdy = (gradQ(ip1,YY)+gradQ(ip2,YY)) / volume_of_two_cells
        Vx = V(jedge,XX)
        Vy = V(jedge,YY)
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(jedge) = abs(aun(jedge))*(Q(ip2)-Q(ip1))*0.5_jprb &
          &          -0.5_jprb*dt*aun(jedge)*(Vx*dQdx+Vy*dQdy)
      end do

      ! non-oscillatory option
      if (limited) then
        call limit_antidiffusive_velocity()
      endif

      ! Compute fluxes from (limited) antidiffusive velocity
      adv(:)   =  0.
      do jedge = 1,nb_edges
        ip1 = edges(jedge,1)
        ip2 = edges(jedge,2)
        flux = aun(jedge)
        adv(ip1) = adv(ip1) + flux
        adv(ip2) = adv(ip2) - flux
      end do

      ! Update the unknowns in vertices
      do jnode=1,nb_nodes
        Q(jnode) = Q(jnode) - adv(jnode)/vol(jnode) * dt
      end do

    end do ! other passes
    end associate

  contains
    
    subroutine compute_Qmax_and_Qmin()
      Qmax(:)  = -1.e10
      Qmin(:)  =  1.e10
      do jedge = 1,geom%nb_edges
        ip1 = geom%edges(jedge,1)
        ip2 = geom%edges(jedge,2)
        Q1 = Q(ip1)
        Q2 = Q(ip2)
        Qmax(ip1) = max( Qmax(ip1), Q1, Q2 )
        Qmin(ip1) = min( Qmin(ip1), Q1, Q2 )
        Qmax(ip2) = max( Qmax(ip2), Q1, Q2 )
        Qmin(ip2) = min( Qmin(ip2), Q1, Q2 )
      end do

      if (Q_is_vector) then
        do jedge = 1,geom%nb_pole_edges
          iedge  = geom%pole_edges(jedge)
          ip1 = geom%edges(iedge,1)
          ip2 = geom%edges(iedge,2)
          Q1 = Q(ip1)
          Q2 = Q(ip2)
          Qmax(ip1) = max( Qmax(ip1), Q1, abs(Q2) )
          Qmin(ip1) = min( Qmin(ip1), Q1,-abs(Q2) )
          Qmax(ip2) = max( Qmax(ip2), Q2, abs(Q1) )
          Qmin(ip2) = min( Qmin(ip2), Q2,-abs(Q1) )
        end do
      end if
    end subroutine compute_Qmax_and_Qmin

    subroutine limit_antidiffusive_velocity
      rhin(:)  =  0.
      rhout(:) =  0.
      do jedge = 1,geom%nb_edges
        ip1 = geom%edges(jedge,1)
        ip2 = geom%edges(jedge,2)
        apos = max(0._jprb,aun(jedge))
        aneg = min(0._jprb,aun(jedge))
        rhin(ip1)  = rhin(ip1)  - aneg
        rhout(ip1) = rhout(ip1) + apos
        rhin(ip2)  = rhin(ip2)  + apos
        rhout(ip2) = rhout(ip2) - aneg
      end do

      do jnode=1,geom%nb_nodes
        cp(jnode) = ( Qmax(jnode)-Q(jnode) )*vol(jnode)/( rhin(jnode) * dt + eps )
        cn(jnode) = ( Q(jnode)-Qmin(jnode) )*vol(jnode)/( rhout(jnode)* dt + eps )
      enddo

      do jedge = 1,geom%nb_edges
        ip1 = geom%edges(jedge,1)
        ip2 = geom%edges(jedge,2)
        if(aun(jedge) > 0._jprb) then
          aun(jedge)=aun(jedge)*min(1._jprb,cp(ip2),cn(ip1))
        else
          aun(jedge)=aun(jedge)*min(1._jprb,cn(ip2),cp(ip1))
        endif
      enddo
    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_gauge

  subroutine compute_gradient(Q,gradQ,Q_is_vector,geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), intent(in)    :: Q(:)
    real(kind=jprb), intent(inout) :: gradQ(:,:)
    logical, intent(in) :: Q_is_vector
    real(kind=jprb), pointer :: S(:,:)
    real(kind=jprb) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2

    S   => vector_field("dual_normals",geom)

    ! derivatives 
    gradQ(:,:) = 0

    do jedge = 1,geom%nb_edges
      ip1 = geom%edges(jedge,1)
      ip2 = geom%edges(jedge,2)
      Sx  = S(jedge,XX)
      Sy  = S(jedge,YY)
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
        iedge = geom%pole_edges(jedge)
        ip1   = geom%edges(iedge,1)
        ip2   = geom%edges(iedge,2)
        Sy   = S(iedge,YY)
        avgQ = ( Q(ip1) + Q(ip2) )*0.5_jprb
        ! correct for wrongly subtracting sy*avgQ in previous loop,
        ! instead of adding, because point at other side of pole
        gradQ(ip2,YY) = gradQ(ip2,YY) + 2*Sy*avgQ 
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
! - setup_shallow_water : To modify dual-volumes to map to the sphere
!                         and create necessary fields
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
  public :: setup_shallow_water 
  public :: propagate_state
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

  subroutine setup_shallow_water(geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb) :: y, G, cos_y, sin_y
    real(kind=jprb), pointer :: coords(:,:), vol(:), hx(:), hy(:), dhxdy_over_G(:), cor(:)
    integer :: jnode
    call create_scalar_field_in_nodes("hx",geom)
    call create_scalar_field_in_nodes("hy",geom)
    call create_scalar_field_in_nodes("dhxdy_over_G",geom)
    call create_scalar_field_in_nodes("coriolis",geom)

    coords => vector_field("coordinates",geom)
    vol    => scalar_field("dual_volumes",geom)
    hx    => scalar_field("hx",geom)
    hy    => scalar_field("hy",geom)
    dhxdy_over_G => scalar_field("dhxdy_over_G",geom)
    cor => scalar_field("coriolis",geom)
    !dir$ ivdep
    do jnode=1,geom%nb_nodes
      y = coords(jnode,YY)
      cos_y = cos(y)
      sin_y = sin(y)
      hx(jnode) = radius*cos_y
      hy(jnode) = radius
      G = hx(jnode)*hy(jnode)
      vol(jnode) = vol(jnode)*G
      dhxdy_over_G(jnode) = - sin_y/(radius*max(eps,cos_y))
      cor(jnode) = f0*sin_y
    enddo

    call create_scalar_field_in_nodes("depth",geom)
    call create_vector_field_in_nodes("momentum",geom)
    call create_vector_field_in_nodes("momentum_forcing",geom)
    call create_scalar_field_in_nodes("depth_backup",geom)
    call create_vector_field_in_nodes("momentum_backup",geom)
    call create_vector_field_in_edges("advective_velocity",geom)
  end subroutine setup_shallow_water


  subroutine propagate_state(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout), target :: geom
    real(kind=jprb) :: tmax, t0, dt_fwd
    character(len=200) :: step_info
    call log_debug( "Propagating state" )
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
      call compute_advective_velocities(dt,geom,"extrapolate")
      call compute_forcing(geom)
    end if
    
    call add_forcing_to_solution(dt,geom)
    call advect_solution(dt,geom)
    call implicit_solve(dt,geom)

    call compute_advective_velocities(dt,geom,"advect")
    call backup_solution(geom)

    geom%fields%time = geom%fields%time+dt
    step = step+1

  end subroutine step_forward




  subroutine set_state_rossby_haurwitz(geom)
    type(DataStructure_type), intent(inout)      :: geom
    real(kind=jprb), dimension(:), pointer   :: D
    real(kind=jprb), dimension(:,:), pointer :: Q, coords
    integer :: jnode, ir
    real(kind=jprb) :: aaa0,zk,om,ph0,g,x,y

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3
    aaa0 = 0.

    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    coords => vector_field("coordinates",geom)

    do jnode=1,geom%nb_nodes
      x=coords(jnode,XX)
      y=coords(jnode,YY)
      if(x == 2._jprb*pi) x=0.
      Q(jnode,XX) =  radius*om*cos(y)+radius*zk*cos(ir*x) *(cos(y))**(ir-1._jprb) &
        &            * (ir*(sin(y))**2-(cos(y))**2)
      Q(jnode,YY) = -radius*zk*ir*(cos(y))**(ir-1._jprb)*sin(y)*sin(ir*x)
      D(jnode) = (ph0+radius**2*fa(y)+radius**2*fb(y)*cos(ir*x) &
        &        +radius**2*fc(y)*cos(2._jprb*ir*x)) / grav
      D(jnode) = max(aaa0,D(jnode))
      Q(jnode,XX) = Q(jnode,XX) * D(jnode)
      Q(jnode,YY) = Q(jnode,YY) * D(jnode)
      if(y == 0.5_jprb*pi) Q(jnode,XX)=0.
      if(y ==-0.5_jprb*pi) Q(jnode,XX)=0.
    end do

    contains 
    ! Helper functions

      real(kind=jprb) function fa(th)
        real(kind=jprb), intent(in) :: th
        fa = om*0.5*(f0+om)*(cos(th))**2 &
          & +0.25*zk**2*(cos(th))**(2*ir)*( (ir+1)*(cos(th))**2 &
          & +(2._jprb*ir**2-ir-2._jprb)-2._jprb*ir**2/(cos(th))**2 )
      end function fa

      real(kind=jprb) function fb(th)
        real(kind=jprb), intent(in) :: th
        fb = (f0+2._jprb*om)*zk/((ir+1._jprb)*(ir+2._jprb))*(cos(th))**ir &
          & *( (ir**2+2._jprb*ir+2._jprb)-((ir+1._jprb)*cos(th))**2 )
      end function fb

      real(kind=jprb) function fc(th)
        real(kind=jprb), intent(in) :: th
        fc = 0.25*zk**2*(cos(th))**(2*ir)*((ir+1._jprb)*(cos(th))**2 -(ir+2._jprb))  
      end function fc

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



  subroutine compute_advective_velocities(dt,geom,option)
    ! this really computes V = G*contravariant_velocity, 
    ! with    G=hx*hy,
    !         physical_velocity = dotproduct( [hx,hy] , contravariant_velocity )
    ! V = (hx*hy) * [u/hx, v/hy] = [u*hy, v*hx]
    ! and hx = r*cos(y)  ,  hy = r
    ! and Q = [ D*u , D*v ]
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    character(len=*), intent(in), optional :: option
    real(kind=jprb) :: Qx, Qy, Q0x, Q0y, Dmod, D0mod, Vx, Vy, Rx, Ry, dVxdx, dVxdy, dVydx, dVydy
    integer :: jnode, jedge, iedge, ip1, ip2
    real(kind=jprb), dimension(:),   pointer :: D, D0, hx, hy, vol
    real(kind=jprb), dimension(:,:), pointer :: Q, Q0, R, Vedges, coords
    real(kind=jprb) :: Vnodes(geom%nb_nodes,2), grad_Vx(geom%nb_nodes,2), grad_Vy(geom%nb_nodes,2)
    coords => vector_field("coordinates",geom)
    Vedges => vector_field("advective_velocity",geom)
    D      => scalar_field("depth",geom)
    D0     => scalar_field("depth_backup",geom)
    Q      => vector_field("momentum",geom)
    Q0     => vector_field("momentum_backup",geom)
    hx     => scalar_field("hx",geom)
    hy     => scalar_field("hy",geom)
    R      => vector_field("momentum_forcing",geom)
    vol    => scalar_field("dual_volumes",geom)
    associate ( nb_nodes => geom%nb_nodes, nb_edges => geom%nb_edges, edges => geom%edges )

    if( option .eq. "advect") then
      do jnode=1,nb_nodes
        Dmod = max(D(jnode), eps)
        Vnodes(jnode,XX)=Q(jnode,XX)/Dmod
        Vnodes(jnode,YY)=Q(jnode,YY)/Dmod
      end do
      call compute_gradient( Vnodes(:,XX), grad_Vx, .True., geom )
      call compute_gradient( Vnodes(:,YY), grad_Vy, .True., geom )

      !dir$ ivdep
      do jnode=1,nb_nodes
        Dmod = max(D(jnode), eps)
        Vx = Vnodes(jnode,XX)
        Vy = Vnodes(jnode,YY)
        Rx = R(jnode,XX)
        Ry = R(jnode,YY)
        dVxdx = grad_Vx(jnode,XX)*hy(jnode)/vol(jnode)
        dVxdy = grad_Vx(jnode,YY)*hx(jnode)/vol(jnode)
        dVydx = grad_Vy(jnode,YY)*hy(jnode)/vol(jnode)    
        dVydy = grad_Vy(jnode,YY)*hx(jnode)/vol(jnode)

        Vnodes(jnode,XX) = ( Vx - 0.5*dt*(Vx*dVxdx+Vy*dVxdy) + 0.5*dt*Rx/Dmod ) * hy(jnode)
        Vnodes(jnode,YY) = ( Vy - 0.5*dt*(Vx*dVydx+Vy*dVydy) + 0.5*dt*Ry/Dmod ) * hx(jnode)
      enddo
    else if( option .eq. "extrapolate") then
      !dir$ ivdep
      do jnode=1,nb_nodes
        Qx    = Q(jnode,XX)
        Qy    = Q(jnode,YY)
        Dmod  = max( eps, D(jnode) )
        Q0x   = Q0(jnode,XX)
        Q0y   = Q0(jnode,YY)
        D0mod = max( eps, D0(jnode) )
        Vnodes(jnode,XX) = ( 1.5_jprb*Qx/Dmod - 0.5_jprb*Q0x/D0mod ) * hy(jnode)
        Vnodes(jnode,YY) = ( 1.5_jprb*Qy/Dmod - 0.5_jprb*Q0y/D0mod ) * hx(jnode)
      end do
    end if

    do jedge=1,nb_edges
      ip1 = edges(jedge,1)
      ip2 = edges(jedge,2)
      Vedges(jedge,XX) = (Vnodes(ip1,XX)+Vnodes(ip2,XX))*0.5_jprb
      Vedges(jedge,YY) = (Vnodes(ip1,YY)+Vnodes(ip2,YY))*0.5_jprb
    enddo
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,geom%nb_pole_edges
      iedge = geom%pole_edges(jedge)
      Vedges(iedge,YY) = 0.
    enddo
    end associate
  end subroutine compute_advective_velocities

  subroutine compute_forcing(geom)
    type(DataStructure_type), intent(inout) :: geom
    integer :: jnode
    real(kind=jprb) :: Qx, Qy, Dmod
    real(kind=jprb), dimension(:),   pointer :: D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprb), dimension(:,:), pointer :: Q, R, coords
    real(kind=jprb) :: grad_D(geom%nb_nodes, 2)
    coords => vector_field("coordinates",geom)
    vol => scalar_field("dual_volumes",geom)
    hx => scalar_field("hx",geom)
    hy => scalar_field("hy",geom)
    dhxdy_over_G => scalar_field("dhxdy_over_G",geom)
    cor => scalar_field("coriolis",geom)

    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)
    call compute_gradient( D, grad_D, .False., geom )
    
    !dir$ ivdep
    do jnode=1,geom%nb_nodes
      Qx    = Q(jnode,XX)
      Qy    = Q(jnode,YY)
      Dmod  = max( eps, D(jnode) )
      R(jnode,XX) = -grav*D(jnode)*grad_D(jnode,XX)*hy(jnode)/vol(jnode) &
        &           + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dmod
      R(jnode,YY) = -grav*D(jnode)*grad_D(jnode,YY)*hx(jnode)/vol(jnode) &
        &           - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dmod
    end do

  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprb), dimension(:,:), pointer :: Q, R
    integer :: jnode, nb_nodes
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)
    !dir$ ivdep
    do jnode=1,geom%nb_nodes
      Q(jnode,XX) = Q(jnode,XX) + 0.5_jprb*dt*R(jnode,XX)
      Q(jnode,YY) = Q(jnode,YY) + 0.5_jprb*dt*R(jnode,YY)
    end do
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,geom)
    real(kind=jprb), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    integer :: jnode, m
    real(kind=jprb) :: Qx, Qy, Rx, Ry, Dmod
    real(kind=jprb) :: Qx_adv, Qy_adv, Rx_exp, Ry_exp

    real(kind=jprb), dimension(:),   pointer :: D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprb), dimension(:,:), pointer :: Q, R, coords
    real(kind=jprb) :: grad_D(geom%nb_nodes, 2)

    coords => vector_field("coordinates",geom)
    vol => scalar_field("dual_volumes",geom)
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)
    hx => scalar_field("hx",geom)
    hy => scalar_field("hy",geom)
    dhxdy_over_G => scalar_field("dhxdy_over_G",geom)
    cor => scalar_field("coriolis",geom)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    call compute_gradient( D, grad_D, .False., geom )
    
    !dir$ ivdep
    do jnode=1,geom%nb_nodes
      Qx    = Q(jnode,XX)
      Qy    = Q(jnode,YY)
      Dmod  = max( eps, D(jnode) )

      Rx_exp = -grav*D(jnode)*grad_D(jnode,XX)*hy(jnode)/vol(jnode)
      Ry_exp = -grav*D(jnode)*grad_D(jnode,YY)*hx(jnode)/vol(jnode)

      Qx_adv = Qx
      Qy_adv = Qy

      do m=1,3 ! Three iterations at most is enough to converge
        Rx = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dmod
        Ry = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dmod
        Qx = Qx_adv + 0.5_jprb*dt*Rx
        Qy = Qy_adv + 0.5_jprb*dt*Ry
      end do
      Q(jnode,XX) = Qx
      Q(jnode,YY) = Qy
      R(jnode,XX) = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dmod
      R(jnode,YY) = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dmod
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

    ! remove noise from periodic boundary... This will dissapear in MPI implementation
    do jnode=1,geom%nb_ghost_nodes
      ip1 = geom%ghost_nodes(jnode,1)
      ip2 = geom%ghost_nodes(jnode,2)
      D(ip2)    = D(ip1)
      Q(ip2,XX) = Q(ip1,XX)
      Q(ip2,YY) = Q(ip1,YY)
    end do
  end subroutine advect_solution

end module shallow_water_module