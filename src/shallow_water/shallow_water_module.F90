! =====================================================================
! mpdata_module
! This module contains strictly algorithmic subroutines
! - mpdata_gauge       : Non-oscillatory variable sign MPDATA advection
! - compute_gradient   : Compute gradient of a scalar array
! =====================================================================
#include "common/assertions.h"
module mpdata_module
  use common_module
  use datastruct_module

  implicit none
  private
  public :: mpdata_gauge
  public :: compute_gradient

  real(kind=jprw), parameter :: eps  = 1.e-6

contains

  subroutine mpdata_gauge(dt,Q,V,order,limited,Q_is_vector_component,geom)
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprw), intent(inout) :: Q(:)
    real(kind=jprw), intent(in) :: V(:,:)
    logical, intent(in) :: Q_is_vector_component, limited
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, ip1,ip2
    real(kind=jprw) :: sx, sy, flux, volume_of_two_cells, dQdx, dQdy, Vx, Vy, apos, aneg, Q1, Q2
    real(kind=jprw) :: Qmin(geom%nb_nodes)
    real(kind=jprw) :: Qmax(geom%nb_nodes)
    real(kind=jprw) :: rhin
    real(kind=jprw) :: rhout
    real(kind=jprw) :: cp(geom%nb_nodes)
    real(kind=jprw) :: cn(geom%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(geom%nb_edges)
    real(kind=jprw) :: fluxv(geom%nb_edges)
    real(kind=jprw) :: gradQ(geom%nb_nodes,2)
    real(kind=jprw), pointer :: vol(:), S(:,:), pole_bc(:)

    vol     => scalar_field("dual_volumes",geom)
    S       => vector_field("dual_normals",geom)
    pole_bc => scalar_field("pole_bc",geom)


    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limited .and. (order .ge. 2) ) then
      if (.not. Q_is_vector_component) call compute_Qmax_and_Qmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,Sx,Sy,Vx,Vy,ip1,ip2,apos,aneg)
    do jedge = 1,geom%nb_edges
      Sx = S(jedge,XX)
      Sy = S(jedge,YY)
      Vx = V(jedge,XX)
      Vy = V(jedge,YY)
      aun(jedge) = Vx*Sx + Vy*Sy
      ip1 = geom%edges(jedge,1)
      ip2 = geom%edges(jedge,2)
      apos = max(0._jprw,aun(jedge))
      aneg = min(0._jprw,aun(jedge))
      fluxv(jedge) = Q(ip1)*apos + Q(ip2)*aneg
    enddo
    !$OMP END PARALLEL DO

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,geom%nb_nodes
      adv = 0.0
      if(geom%nb_neighbours(jnode) > 1) then
        !dir$ unroll 2
        do jedge = 1,geom%nb_neighbours(jnode)
          iedge = geom%my_edges(jedge,jnode)
          adv = adv + geom%sign(jedge,jnode)*fluxv(iedge)
        enddo
      endif
     ! Update the unknowns in vertices
     Q(jnode) = Q(jnode) - adv/vol(jnode) * dt
    enddo
    !$OMP END PARALLEL DO

    call synchronise(Q,geom) ! Qmax and Qmin could be synced here

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector_component, geom)

      call synchronise(gradQ,geom)

      ! Compute antidiffusive normal velocity in faces

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dQdx,dQdy,Vx,Vy)
      do jedge = 1,geom%nb_edges
        ip1 = geom%edges(jedge,1)
        ip2 = geom%edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(ip1) + vol(ip2)
        dQdx = (gradQ(ip1,XX)+gradQ(ip2,XX)) / volume_of_two_cells
        dQdy = (gradQ(ip1,YY)+gradQ(ip2,YY)) / volume_of_two_cells
        Vx = V(jedge,XX)
        Vy = V(jedge,YY)
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(jedge) = abs(aun(jedge))*(Q(ip2)-Q(ip1))*0.5_jprw &
          &          -0.5_jprw*dt*aun(jedge)*(Vx*dQdx+Vy*dQdy)
      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limited) then
        if (Q_is_vector_component) call compute_Qmax_and_Qmin()
        call limit_antidiffusive_velocity()
      endif

      ! Compute fluxes from (limited) antidiffusive velocity
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,adv,jedge,iedge)
      do jnode=1,geom%nb_nodes
        adv = 0.0
        if(geom%nb_neighbours(jnode) > 1) then
          !dir$ unroll 2
          do jedge = 1,geom%nb_neighbours(jnode)
            iedge = geom%my_edges(jedge,jnode)
            adv = adv + geom%sign(jedge,jnode)*aun(iedge)
          enddo
        endif
        ! Update the unknowns in vertices
          Q(jnode) = Q(jnode) - adv/vol(jnode) * dt
      enddo
      !$OMP END PARALLEL DO

      call synchronise(Q,geom)

    end do ! other passes

  contains
    
    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1, Q2    
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jedge,Q1,Q2)
      do jnode=1,geom%nb_nodes
        Q1 = Q(jnode)
        Qmax(jnode) = Q1
        Qmin(jnode) = Q1
        do jedge = 1,geom%nb_neighbours(jnode)
          Q2 = Q(geom%neighbours(jedge,jnode))
          if (.not. Q_is_vector_component) then
            Qmax(jnode) = max( Qmax(jnode), Q2 )
            Qmin(jnode) = min( Qmin(jnode), Q2 )
          else
            Qmax(jnode) = max( Qmax(jnode), pole_bc(jedge)*Q2 )
            Qmin(jnode) = min( Qmin(jnode), pole_bc(jedge)*Q2 )
          end if
        end do
      end do
      !$OMP END PARALLEL DO

      call synchronise(Qmin,geom)
      call synchronise(Qmax,geom)

    end subroutine compute_Qmax_and_Qmin

    subroutine limit_antidiffusive_velocity

      real(kind=jprw) :: asignp,asignn
      real(kind=jprw) :: limit = 1.  ! 1: second order, 0: first order
      if (.not. Q_is_vector_component) limit = 0.995

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,geom%nb_nodes
        rhin  = 0.
        rhout = 0.
        !dir$ unroll 2
        do jedge = 1,geom%nb_neighbours(jnode)
          iedge = geom%my_edges(jedge,jnode)
          apos = max(0._jprw,aun(iedge))
          aneg = min(0._jprw,aun(iedge))
          asignp = max(0._jprw,geom%sign(jedge,jnode))
          asignn = min(0._jprw,geom%sign(jedge,jnode))
          rhin  = rhin  - asignp*aneg - asignn*apos
          rhout = rhout + asignp*apos + asignn*aneg
        end do
        cp(jnode) = ( Qmax(jnode)-Q(jnode) )*vol(jnode)/( rhin * dt + eps )
        cn(jnode) = ( Q(jnode)-Qmin(jnode) )*vol(jnode)/( rhout* dt + eps )
      end do
      !$OMP END PARALLEL DO

      call synchronise(cp,geom)
      call synchronise(cn,geom)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge = 1,geom%nb_edges
        ip1 = geom%edges(jedge,1)
        ip2 = geom%edges(jedge,2)
        if(aun(jedge) > 0._jprw) then
          aun(jedge)=aun(jedge)*min(limit,cp(ip2),cn(ip1))
        else
          aun(jedge)=aun(jedge)*min(limit,cn(ip2),cp(ip1))
        end if
      end do
      !$OMP END PARALLEL DO

    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_gauge

  subroutine compute_gradient(Q,gradQ,Q_is_vector_component,geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprw), intent(in)    :: Q(:)
    real(kind=jprw), intent(inout) :: gradQ(:,:)
    logical, intent(in) :: Q_is_vector_component
    real(kind=jprw), pointer :: S(:,:)
    real(kind=jprw) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2,jnode
    real(kind=jprw) :: avgQSx(geom%nb_edges)
    real(kind=jprw) :: avgQSy(geom%nb_edges)

    S   => vector_field("dual_normals",geom)

    ! derivatives 

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,geom%nb_edges
      ip1 = geom%edges(jedge,1)
      ip2 = geom%edges(jedge,2)
      Sx  = S(jedge,XX)
      Sy  = S(jedge,YY)
      avgQSx(jedge) = Sx*( Q(ip1) + Q(ip2) )*0.5_jprw
      avgQSy(jedge) = Sy*( Q(ip1) + Q(ip2) )*0.5_jprw
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,jedge,iedge)
    do jnode=1,geom%nb_nodes
      gradQ(jnode,XX) = 0.
      gradQ(jnode,YY) = 0.
      do jedge = 1,geom%nb_neighbours(jnode)
        iedge = geom%my_edges(jedge,jnode)
        gradQ(jnode,XX) = gradQ(jnode,XX)+geom%sign(jedge,jnode)*avgQSx(iedge)
        gradQ(jnode,YY) = gradQ(jnode,YY)+geom%sign(jedge,jnode)*avgQSy(iedge)
      end do
    end do
    !$OMP END PARALLEL DO

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    if (.not. Q_is_vector_component) then
      do jedge = 1,geom%nb_pole_edges
        iedge = geom%pole_edges(jedge)
        ip1   = geom%edges(iedge,1)
        ip2   = geom%edges(iedge,2)
        Sy    = S(iedge,YY)
        avgQ  = ( Q(ip1) + Q(ip2) )*0.5_jprw

        ! correct for wrong Y-derivatives in previous loop,
        gradQ(ip2,YY) = gradQ(ip2,YY) + 2._jprw*Sy*avgQ 

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
  use parallel_module
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
  public :: set_state_zonal_flow
  public :: set_time_step

  real(kind=jprw), parameter :: eps    = 1.e-6
  real(kind=jprw), parameter :: radius = 6371.22e+03
  real(kind=jprw), parameter :: f0     = 1.4584e-04 !coriolis parameter (=2xearth's omega)
  real(kind=jprw), parameter :: grav   = 9.80616
  real(kind=jprw), parameter :: pi     = acos(-1._jprw)


  real(kind=jprw) :: dt_forward = 20.
  integer :: iter = 0

contains
 
  subroutine set_time_step(dt)
    real(kind=jprw), intent(in) :: dt
    dt_forward = dt
  end subroutine set_time_step

  subroutine setup_shallow_water(geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprw) :: y, G, cos_y, sin_y
    real(kind=jprw), pointer :: coords(:,:), vol(:), hx(:), hy(:), dhxdy_over_G(:), pole_bc(:)
    integer :: jnode, jedge, iedge
    call create_scalar_field_in_nodes("hx",geom)
    call create_scalar_field_in_nodes("hy",geom)
    call create_scalar_field_in_nodes("dhxdy_over_G",geom)
    call create_scalar_field_in_edges("pole_bc",geom)


    coords => vector_field("coordinates",geom)
    vol    => scalar_field("dual_volumes",geom)
    hx    => scalar_field("hx",geom)
    hy    => scalar_field("hy",geom)
    dhxdy_over_G => scalar_field("dhxdy_over_G",geom)
    pole_bc => scalar_field("pole_bc",geom)
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
    enddo

    pole_bc(:) = 1.
    do jedge=1,geom%nb_pole_edges
      iedge = geom%pole_edges(jedge)
      pole_bc(iedge) = -1.
    end do

    call create_scalar_field_in_nodes("coriolis",geom)
    call create_scalar_field_in_nodes("depth",geom)
    call create_vector_field_in_nodes("momentum",geom)
    call create_vector_field_in_nodes("momentum_forcing",geom)
    call create_scalar_field_in_nodes("depth_backup",geom)
    call create_vector_field_in_nodes("momentum_backup",geom)
    call create_vector_field_in_edges("advective_velocity",geom)

    call create_scalar_field_in_nodes("Dmax_1",geom)
    call create_scalar_field_in_nodes("Dmin_1",geom)
    call create_scalar_field_in_nodes("Dmax_2",geom)
    call create_scalar_field_in_nodes("Dmin_2",geom)
    call create_scalar_field_in_nodes("Dmax_tot",geom)
    call create_scalar_field_in_nodes("Dmin_tot",geom)

  end subroutine setup_shallow_water

  subroutine propagate_state(dt,geom)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout), target :: geom
    real(kind=jprw) :: tend, t0, dt_fwd, tstart
    character(len=200) :: step_info
    real(kind=jprw), pointer :: D(:)

    D => scalar_field("depth",geom)
    tstart   = geom%fields%time
    tend     = geom%fields%time+dt
    
    call log_info(" ")
    write(log_str,'(A,I10,A,I10)') "Propagating from time ",int(tstart)," to time ",int(tend); call log_info()
    do while (geom%fields%time < tend)
      t0 = geom%fields%time
      dt_fwd = min( dt_forward, tend-t0 )
      call step_forward(iter,dt_fwd,geom)

      if( log_level <= LOG_LEVEL_INFO ) then
        if( myproc .eq. 0 ) then
          call progress_bar(geom%fields%time,tstart,tend)
        end if
      else
        write(step_info,'(A6,I8,A12,F9.1,A12,F8.1,A12,E20.15)') "step = ",iter, &
          & "  time = ",geom%fields%time, &
          & "  dt = ",dt_fwd, "Norm ",L2norm(D)
        CALL LOG_INFO( STEP_INFO )
      end if
    end do

  end subroutine propagate_state




  subroutine step_forward(step,dt,geom)
    integer, intent(inout) :: step
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom

    if (step == 0) then ! Pre-compute forcing
      call compute_forcing(geom)

      call backup_solution(geom)
      call compute_advective_velocities(dt,geom,"extrapolate")
    end if
    
    call add_forcing_to_solution(dt,geom)
    call advect_solution(dt,geom)

    call implicit_solve(dt,geom)

    call compute_advective_velocities(dt,geom,"extrapolate")

    call backup_solution(geom)

    geom%fields%time = geom%fields%time+dt
    step = step+1

  end subroutine step_forward




  subroutine set_state_rossby_haurwitz(geom)
    type(DataStructure_type), intent(inout)      :: geom
    real(kind=jprw), dimension(:), pointer   :: D, cor
    real(kind=jprw), dimension(:,:), pointer :: Q, coords
    integer :: jnode, ir
    real(kind=jprw) :: aaa0,zk,om,ph0,g,x,y, sin_y, cos_y

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3
    aaa0 = 0.

    coords => vector_field("coordinates",geom)
    cor => scalar_field("coriolis",geom)
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y,sin_y,cos_y)
    do jnode=1,geom%nb_nodes
      x=coords(jnode,XX)
      y=coords(jnode,YY)
      sin_y = sin(y)
      cos_y = cos(y)
      cor(jnode) = f0*sin_y
      if(x == 2._jprw*pi) x=0.
      Q(jnode,XX) =  radius*om*cos_y+radius*zk*cos(ir*x) * cos_y**(ir-1._jprw) &
        &            * (ir*(sin_y)**2-cos_y**2)
      Q(jnode,YY) = -radius*zk*ir*cos_y**(ir-1._jprw)*sin_y*sin(ir*x)
      D(jnode) = (ph0+radius**2*fa(y)+radius**2*fb(y)*cos(ir*x) &
        &        +radius**2*fc(y)*cos(2._jprw*ir*x)) / grav
      D(jnode) = max(aaa0,D(jnode))
      Q(jnode,XX) = Q(jnode,XX) * D(jnode)
      Q(jnode,YY) = Q(jnode,YY) * D(jnode)
      if(y == 0.5_jprw*pi) Q(jnode,XX)=0.
      if(y ==-0.5_jprw*pi) Q(jnode,XX)=0.
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


  subroutine set_state_zonal_flow(geom)
    type(DataStructure_type), intent(inout)      :: geom
    real(kind=jprw), dimension(:), pointer   :: D, cor
    real(kind=jprw), dimension(:,:), pointer :: Q, coords
    integer :: jnode, ir
    real(kind=jprw) :: x,y
    real(kind=jprw), parameter :: USCAL = 20.
    real(kind=jprw), parameter :: H00 = grav * 8e3
    real(kind=jprw), parameter :: pvel = USCAL/radius
    real(kind=jprw), parameter :: beta = 0.


    coords => vector_field("coordinates",geom)
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    cor => scalar_field("coriolis",geom)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
    do jnode=1,geom%nb_nodes
      x=coords(jnode,XX)
      y=coords(jnode,YY)
      cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
      D(jnode)     = (H00-radius**2*(f0+pvel)*0.5*pvel*(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2)
      D(jnode)     = max(0., D(jnode)/grav)
      Q(jnode,XX)  =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y) * D(jnode)
      Q(jnode,YY)  = -pvel*sin(x)*sin(beta)*radius * D(jnode)
    end do
    !$OMP END PARALLEL DO

  end subroutine set_state_zonal_flow

  subroutine backup_solution(geom)
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprw), dimension(:),   pointer :: D, D0
    real(kind=jprw), dimension(:,:), pointer :: Q, Q0
    integer :: jnode
    
    D  => scalar_field("depth",geom)
    D0 => scalar_field("depth_backup",geom)
    Q  => vector_field("momentum",geom)
    Q0 => vector_field("momentum_backup",geom)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,geom%nb_nodes
      D0(jnode)   = D(jnode)
      Q0(jnode,:) = Q(jnode,:)
    end do
    !$OMP END PARALLEL DO
  end subroutine backup_solution



  subroutine compute_advective_velocities(dt,geom,option)
    ! this really computes V = G*contravariant_velocity, 
    ! with    G=hx*hy,
    !         physical_velocity = dotproduct( [hx,hy] , contravariant_velocity )
    ! V = (hx*hy) * [u/hx, v/hy] = [u*hy, v*hx]
    ! and hx = r*cos(y)  ,  hy = r
    ! and Q = [ D*u , D*v ]
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    character(len=*), intent(in), optional :: option
    real(kind=jprw) :: Qx, Qy, Q0x, Q0y, Dmod, D0mod, Vx, Vy, Rx, Ry, dVxdx, dVxdy, dVydx, dVydy
    integer :: jnode, jedge, iedge, ip1, ip2
    real(kind=jprw), dimension(:),   pointer :: D, D0, hx, hy, vol
    real(kind=jprw), dimension(:,:), pointer :: Q, Q0, R, Vedges, coords
    real(kind=jprw) :: Vnodes(geom%nb_nodes,2), grad_Vx(geom%nb_nodes,2), grad_Vy(geom%nb_nodes,2)

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

    if( option .eq. "advect") then
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Dmod)
      do jnode=1,geom%nb_nodes
        Dmod = max(D(jnode), eps)
        Rx = R(jnode,XX)
        Ry = R(jnode,YY)
        Vnodes(jnode,XX)=(Q(jnode,XX)+0.5*dt*Rx)/Dmod 
        Vnodes(jnode,YY)=(Q(jnode,YY)+0.5*dt*Ry)/Dmod
      end do
      !$OMP END PARALLEL DO
      call compute_gradient( Vnodes(:,XX), grad_Vx, .True., geom )
      call compute_gradient( Vnodes(:,YY), grad_Vy, .True., geom )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Dmod,Vx,Vy,Rx,Ry,dVxdx,dVxdy,dVydx,dVydy)
      do jnode=1,geom%nb_nodes
        Dmod = max(D(jnode), eps)
        Vx = Vnodes(jnode,XX)
        Vy = Vnodes(jnode,YY)
        Rx = R(jnode,XX)
        Ry = R(jnode,YY)
        dVxdx = grad_Vx(jnode,XX)*hy(jnode)/vol(jnode)
        dVxdy = grad_Vx(jnode,YY)*hx(jnode)/vol(jnode)
        dVydx = grad_Vy(jnode,YY)*hy(jnode)/vol(jnode)    
        dVydy = grad_Vy(jnode,YY)*hx(jnode)/vol(jnode)

        Vnodes(jnode,XX) = ( Vx - 0.5*dt*(Vx*dVxdx+Vy*dVxdy)) * hy(jnode)
        Vnodes(jnode,YY) = ( Vy - 0.5*dt*(Vx*dVydx+Vy*dVydy)) * hx(jnode)
      enddo
      !$OMP END PARALLEL DO
      call synchronise( Vnodes, geom )
    else if( option .eq. "extrapolate") then
      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Qx,Qy,Dmod,Q0x,Q0y,D0mod)
      do jnode=1,geom%nb_nodes
        Qx    = Q(jnode,XX)
        Qy    = Q(jnode,YY)
        Dmod  = max( eps, D(jnode) )
        Q0x   = Q0(jnode,XX)
        Q0y   = Q0(jnode,YY)
        D0mod = max( eps, D0(jnode) )
        Vnodes(jnode,XX) = ( 1.5_jprw*Qx/Dmod - 0.5_jprw*Q0x/D0mod ) * hy(jnode)
        Vnodes(jnode,YY) = ( 1.5_jprw*Qy/Dmod - 0.5_jprw*Q0y/D0mod ) * hx(jnode)
      end do
      !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
    do jedge=1,geom%nb_edges
      ip1 = geom%edges(jedge,1)
      ip2 = geom%edges(jedge,2)
      Vedges(jedge,XX) = (Vnodes(ip1,XX)+Vnodes(ip2,XX))*0.5_jprw
      Vedges(jedge,YY) = (Vnodes(ip1,YY)+Vnodes(ip2,YY))*0.5_jprw
    enddo
    !$OMP END PARALLEL DO
    
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,geom%nb_pole_edges
      iedge = geom%pole_edges(jedge)
      Vedges(iedge,YY) = 0.
    enddo

  end subroutine compute_advective_velocities

  subroutine compute_forcing(geom)
    type(DataStructure_type), intent(inout) :: geom
    integer :: jnode
    real(kind=jprw) :: Qx, Qy, Dmod
    real(kind=jprw), dimension(:),   pointer :: D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:), pointer :: Q, R, coords
    real(kind=jprw) :: grad_D(geom%nb_nodes, 2)
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
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Qx,Qy,Dmod)
    do jnode=1,geom%nb_nodes
      Qx    = Q(jnode,XX)
      Qy    = Q(jnode,YY)
      Dmod  = max( eps, D(jnode) )
      R(jnode,XX) = -grav*D(jnode)*grad_D(jnode,XX)*hy(jnode)/vol(jnode) &
        &           + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dmod
      R(jnode,YY) = -grav*D(jnode)*grad_D(jnode,YY)*hx(jnode)/vol(jnode) &
        &           - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dmod
    end do
    !$OMP END PARALLEL DO

    call synchronise(R,geom)
  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,geom)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprw), dimension(:,:), pointer :: Q, R
    integer :: jnode, nb_nodes
    Q => vector_field("momentum",geom)
    R => vector_field("momentum_forcing",geom)
    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,geom%nb_nodes
      Q(jnode,XX) = Q(jnode,XX) + 0.5_jprw*dt*R(jnode,XX)
      Q(jnode,YY) = Q(jnode,YY) + 0.5_jprw*dt*R(jnode,YY)
    end do
    !$OMP END PARALLEL DO
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,geom)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    integer :: jnode, m
    real(kind=jprw) :: Qx, Qy, Rx, Ry, Dmod
    real(kind=jprw) :: Qx_adv, Qy_adv, Rx_exp, Ry_exp

    real(kind=jprw), dimension(:),   pointer :: D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:), pointer :: Q, R, coords
    real(kind=jprw) :: grad_D(geom%nb_nodes, 2)

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
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Qx,Qy,Dmod,Rx_exp,Ry_exp,Qx_adv,Qy_adv,m,Rx,Ry)
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
        Qx = Qx_adv + 0.5_jprw*dt*Rx
        Qy = Qy_adv + 0.5_jprw*dt*Ry
      end do
      Q(jnode,XX) = Qx
      Q(jnode,YY) = Qy
      R(jnode,XX) = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dmod
      R(jnode,YY) = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dmod
    end do
    !$OMP END PARALLEL DO
    call synchronise(Q,geom)
    call synchronise(R,geom)

  end subroutine implicit_solve



  subroutine advect_solution(dt,geom)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: geom
    real(kind=jprw), dimension(:),   pointer :: D
    real(kind=jprw), dimension(:,:), pointer :: Q, V
    integer :: jnode, ip1, ip2
    D => scalar_field("depth",geom)
    Q => vector_field("momentum",geom)
    V => vector_field("advective_velocity",geom)

    !    mpdata_gauge( time, variable, velocity, order, limit,  is_vector, geom )
    call mpdata_gauge( dt,   D,        V,        2,     .True., .False.,   geom )
    call mpdata_gauge( dt,   Q(:,XX),  V,        2,     .True., .True. ,   geom )
    call mpdata_gauge( dt,   Q(:,YY),  V,        2,     .True., .True. ,   geom )
  end subroutine advect_solution

end module shallow_water_module
