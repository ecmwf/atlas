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
  public :: mpdata_D
  public :: mpdata_gauge_U
  public :: mpdata_Q
  public :: compute_gradient
  public :: compute_gradient_tensor

 integer, parameter, public :: MPDATA_STANDARD=1, MPDATA_GAUGE=2

  real(kind=jprw), parameter :: eps = 1.e-6

  integer, parameter, public :: probe = 22186
  real(kind=jprw), parameter :: Dphys_min = -10.
  real(kind=jprw), parameter :: Dphys_max = 20000


contains

  subroutine mpdata_D(mpdata_scheme,dt,D,V,VDS,order,limit,dstruct)
    ! For mpdata standard scheme, in the case for positive D,
    !  grad( abs(D) ) == grad( D )
    integer, intent(in) :: mpdata_scheme
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: D(:)
    real(kind=jprw), intent(in) :: V(:,:)
    real(kind=jprw), intent(in) :: limit
    integer, intent(in) :: order
    real(kind=jprw), intent(out) :: VDS(dstruct%nb_edges)

    integer :: jnode, jedge, iedge, jpass, ip1,ip2
    real(kind=jprw) :: sx, sy, volume_of_two_cells, dDdx, dDdy, Vx, Vy, apos, aneg, Dtmp, D_abs
    real(kind=jprw) :: Dmin(dstruct%nb_nodes)
    real(kind=jprw) :: Dmax(dstruct%nb_nodes)
    real(kind=jprw) :: rhin
    real(kind=jprw) :: rhout
    real(kind=jprw) :: cp(dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(dstruct%nb_edges)
    real(kind=jprw) :: gradD(2,dstruct%nb_nodes)
    real(kind=jprw) :: sum_Sabs(2,dstruct%nb_nodes)
    real(kind=jprw) :: sum_Dbar(2,dstruct%nb_nodes)
    real(kind=jprw), pointer :: vol(:), S(:,:), Wind(:,:)
    real(kind=jprw) :: Dbar(2) ! Dabs_bar == Dbar since D > 0 always

    Wind    => vector_field_2d("velocity",dstruct)

    vol     => scalar_field_2d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)

    VDS(:) = 0.

    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limit > 0._jprw .and. (order >= 2) ) then
      Dmax(:) = -1e10
      Dmin(:) =  1e10
      call compute_Dmax_and_Dmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,Sx,Sy,Vx,Vy,ip1,ip2,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

      Sx = S(XX,jedge)
      Sy = S(YY,jedge)
      Vx = V(XX,jedge)
      Vy = V(YY,jedge)

      aun(jedge) = Vx*Sx + Vy*Sy

      apos = max(0._jprw,aun(jedge))
      aneg = min(0._jprw,aun(jedge))
      VDS(jedge) = D(ip1)*apos + D(ip2)*aneg
    enddo
    !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          adv = adv + dstruct%sign(jedge,jnode)*VDS(iedge)
        enddo
      endif
      ! Update the unknowns in vertices
      Dtmp = D(jnode) - adv/vol(jnode) * dt

      if (Dtmp < Dphys_min .or. Dtmp > Dphys_max ) then
        write(log_str,*) "D_pass0 ", jnode, D(jnode); call log_error()
        write(log_str,*) "D_pass1 ", jnode, Dtmp; call log_error()
        write(log_str,*) "adv     ", jnode, adv; call log_error()
        write(log_str,*) "velocity", jnode, Wind(:,jnode); call log_error()
        call abort
      end if
      if ( jnode .eq. probe ) then
        write(log_str,*) "D_pass0", jnode, D(jnode); call log_debug()
        write(log_str,*) "D_pass1", jnode, Dtmp; call log_debug()
        write(log_str,*) "adv     ", jnode, adv; call log_debug()
        write(log_str,*) "velocity", jnode, Wind(:,jnode); call log_debug()
      end if

      D(jnode) = max(eps, Dtmp)

    enddo

    !$OMP END PARALLEL DO

    call halo_exchange(D,dstruct) ! Dmax and Dmin could be synced here

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------

    do jpass=2,order

      ! Compute derivatives for mpdata
      select case (mpdata_scheme)
        case (MPDATA_STANDARD)
          call compute_gradient_abs(D, gradD, dstruct)
        case (MPDATA_GAUGE)
          call compute_gradient(D, gradD, dstruct)
      end select

      call halo_exchange(gradD,dstruct)

      ! Compute antidiffusive normal velocity in faces

      if (mpdata_scheme == MPDATA_STANDARD) then
        sum_Dbar(:,:) = 0._jprw
        sum_Sabs(:,:) = 0._jprw
        do jedge=1, dstruct%nb_edges
          ip1 = dstruct%edges(jedge,1)
          ip2 = dstruct%edges(jedge,2)
          Sx = abs(S(XX,jedge))
          Sy = abs(S(YY,jedge))
          D_abs = ( abs(D(ip1)) + abs(D(ip2)) ) * 0.5_jprw
          sum_Dbar(XX,ip1) = sum_Dbar(XX,ip1) + Sx * D_abs
          sum_Dbar(XX,ip2) = sum_Dbar(XX,ip2) + Sx * D_abs
          sum_Dbar(YY,ip1) = sum_Dbar(YY,ip1) + Sy * D_abs
          sum_Dbar(YY,ip2) = sum_Dbar(YY,ip2) + Sy * D_abs
          sum_Sabs(XX,ip1) = sum_Sabs(XX,ip1) + Sx
          sum_Sabs(XX,ip2) = sum_Sabs(XX,ip2) + Sx
          sum_Sabs(YY,ip1) = sum_Sabs(YY,ip1) + Sy
          sum_Sabs(YY,ip2) = sum_Sabs(YY,ip2) + Sy
        end do
      end if

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dDdx,dDdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(ip1) + vol(ip2)
        dDdx = (gradD(XX,ip1)+gradD(XX,ip2)) / volume_of_two_cells
        dDdy = (gradD(YY,ip1)+gradD(YY,ip2)) / volume_of_two_cells
        Vx = V(XX,jedge)
        Vy = V(YY,jedge)

        select case (mpdata_scheme)
          case (MPDATA_STANDARD)
            Dbar(XX) = (sum_Dbar(XX,ip1) + sum_Dbar(XX,ip2) + eps) / ( sum_Sabs(XX,ip1) + sum_Sabs(XX,ip2) )
            Dbar(YY) = (sum_Dbar(YY,ip1) + sum_Dbar(YY,ip2) + eps) / ( sum_Sabs(YY,ip1) + sum_Sabs(YY,ip2) )
            aun(jedge) = abs(aun(jedge))*( abs(D(ip2))-abs(D(ip1)) )/(abs(D(ip2))+abs(D(ip1))+eps) &
              &          -0.5_jprw*dt*aun(jedge)*(Vx*dDdx/Dbar(XX)+Vy*dDdy/Dbar(YY))
          case (MPDATA_GAUGE)
            ! variable sign option with asymptotic analysis, (mpdata gauge)
            aun(jedge) = abs(aun(jedge))*(D(ip2)-D(ip1))*0.5_jprw &
              &          -0.5_jprw*dt*aun(jedge)*(Vx*dDdx+Vy*dDdy)
        end select

      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limit > 0._jprw) then
        !if (mpdata_scheme == MPDATA_STANDARD) call compute_Dmax_and_Dmin()
        call limit_antidiffusive_velocity()
      endif

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge)
      do jedge = 1,dstruct%nb_edges
        VDS(jedge) = VDS(jedge) + aun(jedge)
      end do
      !$OMP END PARALLEL DO

      ! Compute fluxes from (limited) antidiffusive velocity
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
      do jnode=1,dstruct%nb_nodes
        adv = 0.0
        if(dstruct%nb_neighbours(jnode) > 1) then
          do jedge = 1,dstruct%nb_neighbours(jnode)
            iedge = dstruct%my_edges(jedge,jnode)
            adv = adv + dstruct%sign(jedge,jnode)*aun(iedge)
          enddo
        endif
        ! Update the unknowns in vertices
        !D(jnode) = max( D(jnode) - adv/vol(jnode) * dt, eps )
        Dtmp = D(jnode) - adv/vol(jnode) * dt
        if ( Dtmp < Dphys_min .or. Dtmp > Dphys_max ) then
          write(0,*) "ERROR in node",jnode
          write(0,*) "D_pass2 = ", Dtmp
          write(0,*) "gradD = ", gradD(:,jnode)
          call abort
        end if
        if ( jnode == probe ) then
          write(log_str,*) "D_pass2", jnode, Dtmp; call log_debug()
          write(log_str,*) "gradD = ", gradD(:,jnode); call log_debug()
        end if

        D(jnode) = Dtmp

      enddo
      !$OMP END PARALLEL DO
      call halo_exchange(D,dstruct)


    end do ! other passes

  contains

    subroutine compute_Dmax_and_Dmin( )
      real(kind=jprw) :: D1, D2
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,D1,D2)
      do jnode=1,dstruct%nb_nodes
        D1 = D(jnode)
        Dmax(jnode) = max( Dmax(jnode), D1 )
        Dmin(jnode) = min( Dmin(jnode), D1 )
        do jedge = 1,dstruct%nb_neighbours(jnode)
          D2 = D(dstruct%neighbours(jedge,jnode))
          iedge = dstruct%my_edges(jedge,jnode)
          Dmax(jnode) = max( Dmax(jnode), D2 )
          Dmin(jnode) = min( Dmin(jnode), D2 )
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange(Dmin,dstruct)
      call halo_exchange(Dmax,dstruct)

    end subroutine compute_Dmax_and_Dmin

    subroutine limit_antidiffusive_velocity

      real(kind=jprw) :: asignp,asignn

      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,dstruct%nb_nodes
        rhin  = 0.
        rhout = 0.
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          apos = max(0._jprw,aun(iedge))
          aneg = min(0._jprw,aun(iedge))
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          rhin  = rhin  - asignp*aneg - asignn*apos
          rhout = rhout + asignp*apos + asignn*aneg
        end do
        cp(jnode) = ( Dmax(jnode)-D(jnode) )*vol(jnode)/( rhin * dt + eps )
        cn(jnode) = ( D(jnode)-Dmin(jnode) )*vol(jnode)/( rhout* dt + eps )
      end do
      !$OMP END PARALLEL DO

      call halo_exchange(cp,dstruct)
      call halo_exchange(cn,dstruct)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        if(aun(jedge) > 0._jprw) then
          aun(jedge)=aun(jedge)*min(limit,cp(ip2),cn(ip1))
        else
          aun(jedge)=aun(jedge)*min(limit,cn(ip2),cp(ip1))
        end if
      end do
      !$OMP END PARALLEL DO

    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_D
  
  
  
  
  subroutine mpdata_gauge_U(dt,U,VDS,DR,D,order,limited,dstruct)
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: U(:,:), D(:)
    real(kind=jprw), intent(in) :: VDS(:), DR(:)
    logical, intent(in) :: limited
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, ip1,ip2
    real(kind=jprw) :: Sx, Sy, Ssqr, volume_of_two_cells, dUdx(2), dUdy(2), Vx, Vy
    real(kind=jprw) :: apos(2), aneg(2), x1, x2, y1, y2
    real(kind=jprw) :: Umin(2,dstruct%nb_nodes)
    real(kind=jprw) :: Umax(2,dstruct%nb_nodes)
    real(kind=jprw) :: rhin(2)
    real(kind=jprw) :: rhout(2)
    real(kind=jprw) :: cp(2,dstruct%nb_nodes)
    real(kind=jprw) :: cn(2,dstruct%nb_nodes)
    real(kind=jprw) :: adv(2)
    real(kind=jprw) :: aun(2,dstruct%nb_edges)
    real(kind=jprw) :: fluxv(2,dstruct%nb_edges)
    real(kind=jprw) :: gradU(4,dstruct%nb_nodes)
    real(kind=jprw) :: VD(2,dstruct%nb_edges)
    real(kind=jprw) :: VDnodes(2,dstruct%nb_nodes)
    real(kind=jprw) :: sum_inv_distance
    real(kind=jprw) :: inv_distance(dstruct%nb_nodes)
    real(kind=jprw) :: Utmp(2,dstruct%nb_nodes)
    real(kind=jprw) :: volD(dstruct%nb_nodes)
    real(kind=jprw), pointer :: vol(:), S(:,:), pole_bc(:), coords(:,:)

    vol     => scalar_field_2d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)
    coords  => vector_field_2d("coordinates",dstruct)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy,Ssqr,x1,x2,y1,y2)
    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx = S(XX,jedge)
      Sy = S(YY,jedge)
      Ssqr =  Sx*Sx + Sy*Sy 
      VD(XX,jedge) = Sx/Ssqr * VDS(jedge)
      VD(YY,jedge) = Sy/Ssqr * VDS(jedge)
      x1 = coords(XX,ip1)
      x2 = coords(XX,ip2)
      y1 = coords(YY,ip1)
      y2 = coords(YY,ip2)
      inv_distance(jedge) = 1._jprw/sqrt( (x2-x1)**2 + (y2-y1)**2 )
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,sum_inv_distance,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      volD(jnode) = vol(jnode)*D(jnode)
      sum_inv_distance = 0._jprw
      VDnodes(:,jnode) = 0._jprw
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        sum_inv_distance = sum_inv_distance + inv_distance(iedge)
        VDnodes(:,jnode) = VDnodes(:,jnode) + VD(:,iedge) * inv_distance(iedge)
      enddo
      VDnodes(:,jnode) = VDnodes(:,jnode) / sum_inv_distance
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      VD(:,jedge) = 0.5_jprw * (VDnodes(:,ip1) + VDnodes(:,ip2) )
    end do
    !$OMP END PARALLEL DO


    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      VD(YY,iedge) = 0.
    enddo


    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limited .and. (order .ge. 2) ) then
      Umax(:,:) = -1e10
      Umin(:,:) =  1e10
      call compute_Umax_and_Umin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      aun(:,jedge) = VDS(jedge)
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      apos(:) = max(0._jprw,aun(:,jedge))
      aneg(:) = min(0._jprw,aun(:,jedge))
      fluxv(:,jedge) = U(:,ip1)*apos(:) + U(:,ip2)*aneg(:)
    enddo
    !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv(:) = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          adv(:) = adv(:) + dstruct%sign(jedge,jnode)*fluxv(:,iedge)
        enddo
      endif
     ! Update the unknowns in vertices
     U(:,jnode) = U(:,jnode) - adv(:)/max(eps,volD(jnode)) * dt

    enddo
    !$OMP END PARALLEL DO

    call halo_exchange(U,dstruct) 

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode = 1,dstruct%nb_nodes
      Utmp(:,jnode) = U(:,jnode)
      U(:,jnode) = U(:,jnode)*DR(jnode)
    enddo
    !$OMP END PARALLEL DO

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient_tensor(U, gradU, dstruct)

      call halo_exchange(gradU,dstruct)

      ! Compute antidiffusive normal velocity in faces

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dUdx,dUdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = max(eps, volD(ip1) + volD(ip2) )
        dUdx(XX) = (gradU(XXDXX,ip1)+gradU(XXDXX,ip2)) / volume_of_two_cells
        dUdx(YY) = (gradU(YYDXX,ip1)+gradU(YYDXX,ip2)) / volume_of_two_cells
        dUdy(XX) = (gradU(XXDYY,ip1)+gradU(XXDYY,ip2)) / volume_of_two_cells
        dUdy(YY) = (gradU(YYDYY,ip1)+gradU(YYDYY,ip2)) / volume_of_two_cells
        Vx = VD(XX,jedge)
        Vy = VD(YY,jedge)
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(:,jedge) = abs(aun(:,jedge))*(U(:,ip2)-U(:,ip1))*0.5_jprw &
          &          -0.5_jprw*dt*aun(:,jedge)*(Vx*dUdx(:)+Vy*dUdy(:))  ! = VDS*dUds
      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limited) then
        !call compute_Umax_and_Umin()
        call limit_antidiffusive_velocity()
      endif


      ! Compute fluxes from (limited) antidiffusive velocity
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
      do jnode=1,dstruct%nb_nodes
        adv(:) = 0.0
        if(dstruct%nb_neighbours(jnode) > 1) then
          do jedge = 1,dstruct%nb_neighbours(jnode)
            iedge = dstruct%my_edges(jedge,jnode)
            adv(:) = adv(:) + dstruct%sign(jedge,jnode)*aun(:,iedge)
          enddo
        endif
        ! Update the unknowns in vertices
        U(:,jnode) = Utmp(:,jnode) - adv(:)/max(eps,volD(jnode)) * dt
        U(:,jnode) = U(:,jnode) * DR(jnode)
      enddo
      !$OMP END PARALLEL DO

      call halo_exchange(U,dstruct)

    end do ! other passes

  contains
    
    subroutine compute_Umax_and_Umin( )
      real(kind=jprw) :: U1(2), U2(2)    
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,U1,U2)
      do jnode=1,dstruct%nb_nodes
        U1(:) = U(:,jnode)
        Umax(:,jnode) = max( Umax(:,jnode), U1(:) )
        Umin(:,jnode) = min( Umin(:,jnode), U1(:) )
        do jedge = 1,dstruct%nb_neighbours(jnode)
          U2(:) = U(:,dstruct%neighbours(jedge,jnode))
          iedge = dstruct%my_edges(jedge,jnode)
          Umax(:,jnode) = max( Umax(:,jnode), pole_bc(iedge)*U2(:) )
          Umin(:,jnode) = min( Umin(:,jnode), pole_bc(iedge)*U2(:) )
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange(Umin,dstruct)
      call halo_exchange(Umax,dstruct)

    end subroutine compute_Umax_and_Umin

    subroutine limit_antidiffusive_velocity

      real(kind=jprw) :: asignp,asignn
      real(kind=jprw) :: limit = 1.  ! 1: second order, 0: first order
      integer :: var

      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,dstruct%nb_nodes
        rhin(:)  = 0.
        rhout(:) = 0.
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          apos(:) = max(0._jprw,aun(:,iedge))
          aneg(:) = min(0._jprw,aun(:,iedge))
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          rhin(:)  = rhin(:)  - asignp*aneg(:) - asignn*apos(:)
          rhout(:) = rhout(:) + asignp*apos(:) + asignn*aneg(:)
        end do
        cp(:,jnode) = ( Umax(:,jnode)-U(:,jnode) )*volD(jnode)/( rhin(:) * dt + eps )
        cn(:,jnode) = ( U(:,jnode)-Umin(:,jnode) )*volD(jnode)/( rhout(:)* dt + eps )
      end do
      !$OMP END PARALLEL DO

      call halo_exchange(cp,dstruct)
      call halo_exchange(cn,dstruct)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,var)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do var=1,2
          if(aun(var,jedge) > 0._jprw) then
            aun(var,jedge)=aun(var,jedge)*min(limit,cp(var,ip2),cn(var,ip1))
          else
            aun(var,jedge)=aun(var,jedge)*min(limit,cn(var,ip2),cp(var,ip1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_gauge_U
  
  
  subroutine mpdata_Q(mpdata_scheme,dt,Q,V,order,limit,dstruct)
    integer, intent(in) :: mpdata_scheme
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: Q(:,:)
    real(kind=jprw), intent(in) :: V(:,:)
    real(kind=jprw), intent(in) :: limit
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, ip1,ip2
    real(kind=jprw) :: Sx, Sy, volume_of_two_cells, dQdx(2), dQdy(2), Vx, Vy, Qtmp(2), Q_abs(2), Qbar(2,2)
    real(kind=jprw) :: apos(2), aneg(2)
    real(kind=jprw) :: Qmin(2,dstruct%nb_nodes)
    real(kind=jprw) :: Qmax(2,dstruct%nb_nodes)
    real(kind=jprw) :: rhin(2)
    real(kind=jprw) :: rhout(2)
    real(kind=jprw) :: cp(2,dstruct%nb_nodes)
    real(kind=jprw) :: cn(2,dstruct%nb_nodes)
    real(kind=jprw) :: adv(2)
    real(kind=jprw) :: aun(2,dstruct%nb_edges)
    real(kind=jprw) :: fluxv(2,dstruct%nb_edges)
    real(kind=jprw) :: gradQ(4,dstruct%nb_nodes)
    real(kind=jprw) :: sum_Sabs(2,dstruct%nb_nodes)
    real(kind=jprw) :: sum_Qbar(2,2,dstruct%nb_nodes)
    real(kind=jprw), pointer :: vol(:), S(:,:), pole_bc(:)

    vol     => scalar_field_2d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)


    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limit > 0. .and. (order >= 2) ) then
      Qmax(:,:) = -1e10
      Qmin(:,:) =  1e10
      call compute_Qmax_and_Qmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy,Vx,Vy,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

      Sx = S(XX,jedge)
      Sy = S(YY,jedge)

      Vx = V(XX,jedge)
      Vy = V(YY,jedge)

      aun(:,jedge) = Vx*Sx + Vy*Sy
      apos(:) = max(0._jprw,aun(:,jedge))
      aneg(:) = min(0._jprw,aun(:,jedge))
      fluxv(:,jedge) = Q(:,ip1)*apos(:) + Q(:,ip2)*aneg(:)
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv(:) = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          adv(:) = adv(:) + dstruct%sign(jedge,jnode)*fluxv(:,iedge)
        end do
      end if
     ! Update the unknowns in vertices
     Q(:,jnode) = Q(:,jnode) - adv(:)/vol(jnode) * dt
    end do
    !$OMP END PARALLEL DO

    do jnode=1,dstruct%nb_nodes
      if ( abs(Q(XX,jnode)) > 1e30 .or. abs(Q(YY,jnode)) > 1e30 ) then
        write(0,*) "Q_pass1", jnode, Q(:,jnode)
        call abort
      end if
      if ( jnode .eq. probe ) then
        write(log_str,*) "Q_pass1", jnode, Q(:,jnode); call log_debug()
      end if
    end do

    call halo_exchange(Q,dstruct)


    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------

    do jpass=2,order

      ! Compute derivatives for mpdata
      select case (mpdata_scheme)
        case (MPDATA_STANDARD)
          call compute_gradient_tensor_abs(Q, gradQ, dstruct)
        case (MPDATA_GAUGE)
          call compute_gradient_tensor(Q, gradQ, dstruct)
      end select

      call halo_exchange(gradQ,dstruct)

      ! Compute antidiffusive normal velocity in faces

      if (mpdata_scheme == MPDATA_STANDARD) then
        sum_Qbar(:,:,:) = 0._jprw
        sum_Sabs(:,:) = 0._jprw
        do jedge=1, dstruct%nb_edges
          ip1 = dstruct%edges(jedge,1)
          ip2 = dstruct%edges(jedge,2)
          Sx = abs(S(XX,jedge))
          Sy = abs(S(YY,jedge))
          Q_abs(:) = ( abs(Q(:,ip1)) + abs(Q(:,ip2)) ) * 0.5_jprw
          sum_Qbar(:,XX,ip1) = sum_Qbar(:,XX,ip1) + Sx * Q_abs(:)
          sum_Qbar(:,YY,ip1) = sum_Qbar(:,YY,ip1) + Sy * Q_abs(:)
          sum_Qbar(:,XX,ip2) = sum_Qbar(:,XX,ip2) + Sx * Q_abs(:)
          sum_Qbar(:,YY,ip2) = sum_Qbar(:,YY,ip2) + Sy * Q_abs(:)
          sum_Sabs(XX,ip1) = sum_Sabs(XX,ip1) + Sx
          sum_Sabs(XX,ip2) = sum_Sabs(XX,ip2) + Sx
          sum_Sabs(YY,ip1) = sum_Sabs(YY,ip1) + Sy
          sum_Sabs(YY,ip2) = sum_Sabs(YY,ip2) + Sy
        end do
      end if

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dQdx,dQdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = max(eps, vol(ip1) + vol(ip2) )
        dQdx(XX) = (gradQ(XXDXX,ip1)+gradQ(XXDXX,ip2)) / volume_of_two_cells
        dQdx(YY) = (gradQ(YYDXX,ip1)+gradQ(YYDXX,ip2)) / volume_of_two_cells
        dQdy(XX) = (gradQ(XXDYY,ip1)+gradQ(XXDYY,ip2)) / volume_of_two_cells
        dQdy(YY) = (gradQ(YYDYY,ip1)+gradQ(YYDYY,ip2)) / volume_of_two_cells
        Vx = V(XX,jedge)
        Vy = V(YY,jedge)

        select case (mpdata_scheme)
          case (MPDATA_STANDARD)
            Qbar(:,XX) = (sum_Qbar(:,XX,ip1) + sum_Qbar(:,XX,ip2) + eps) / ( sum_Sabs(XX,ip1) + sum_Sabs(XX,ip2) )
            Qbar(:,YY) = (sum_Qbar(:,YY,ip1) + sum_Qbar(:,YY,ip2) + eps) / ( sum_Sabs(YY,ip1) + sum_Sabs(YY,ip2) )
            aun(:,jedge) = abs(aun(:,jedge))*( abs(Q(:,ip2))-abs(Q(:,ip1)) )/(abs(Q(:,ip2))+abs(Q(:,ip1))+eps) &
              &          -0.5_jprw*dt*aun(:,jedge)*(Vx*dQdx(:)/Qbar(:,XX)+Vy*dQdy(:)/Qbar(:,YY))
          case (MPDATA_GAUGE)
            ! variable sign option with asymptotic analysis, (mpdata gauge)
            aun(:,jedge) = abs(aun(:,jedge))*(Q(:,ip2)-Q(:,ip1))*0.5_jprw &
              &          -0.5_jprw*dt*aun(:,jedge)*(Vx*dQdx(:)+Vy*dQdy(:))
        end select

      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limit > 0._jprw) then
        !if (mpdata_scheme == MPDATA_STANDARD) call compute_Qmax_and_Qmin()
        call limit_antidiffusive_velocity()
      endif


      ! Compute fluxes from (limited) antidiffusive velocity
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
      do jnode=1,dstruct%nb_nodes
        adv(:) = 0.0
        if(dstruct%nb_neighbours(jnode) > 1) then
          do jedge = 1,dstruct%nb_neighbours(jnode)
            iedge = dstruct%my_edges(jedge,jnode)
            adv(:) = adv(:) + dstruct%sign(jedge,jnode)*aun(:,iedge)
          enddo
        endif
        ! Update the unknowns in vertices
        !Q(:,jnode) = Q(:,jnode) - adv(:)/vol(jnode) * dt
        Qtmp(:) = Q(:,jnode) - adv(:)/vol(jnode) * dt
        if ( abs(Qtmp(XX) ) > 1e10 .or. abs(Qtmp(YY) ) > 1e10 ) then
          write(0,*) "Q_pass2 = ", Qtmp(:)
          write(0,*) "gradQ = ", gradQ(:,jnode)
          call abort
        end if
        if ( jnode .eq. probe ) then
          write(log_str,*) "Q_pass2", jnode, Qtmp; call log_debug()
          write(log_str,*) "gradQ  ", jnode, gradQ(:,jnode); call log_debug()
        end if

        Q(:,jnode) = Qtmp(:)
      enddo
      !$OMP END PARALLEL DO

      call halo_exchange(Q,dstruct)

    end do ! other passes

  contains

    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1(2), Q2(2)
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,Q1,Q2)
      do jnode=1,dstruct%nb_nodes
        Q1(:) = Q(:,jnode)
        Qmax(:,jnode) = max( Qmax(:,jnode), Q1(:) )
        Qmin(:,jnode) = min( Qmin(:,jnode), Q1(:) )
        do jedge = 1,dstruct%nb_neighbours(jnode)
          Q2(:) = Q(:,dstruct%neighbours(jedge,jnode))
          iedge = dstruct%my_edges(jedge,jnode)
          Qmax(:,jnode) = max( Qmax(:,jnode), pole_bc(iedge)*Q2(:) )
          Qmin(:,jnode) = min( Qmin(:,jnode), pole_bc(iedge)*Q2(:) )
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange(Qmin,dstruct)
      call halo_exchange(Qmax,dstruct)

    end subroutine compute_Qmax_and_Qmin

    subroutine limit_antidiffusive_velocity

      real(kind=jprw) :: asignp,asignn
      integer :: var
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,dstruct%nb_nodes
        rhin(:)  = 0.
        rhout(:) = 0.
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          apos(:) = max(0._jprw,aun(:,iedge))
          aneg(:) = min(0._jprw,aun(:,iedge))
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          rhin(:)  = rhin(:)  - asignp*aneg(:) - asignn*apos(:)
          rhout(:) = rhout(:) + asignp*apos(:) + asignn*aneg(:)
        end do
        cp(:,jnode) = ( Qmax(:,jnode)-Q(:,jnode) )*vol(jnode)/( rhin(:) * dt + eps )
        cn(:,jnode) = ( Q(:,jnode)-Qmin(:,jnode) )*vol(jnode)/( rhout(:)* dt + eps )
      end do
      !$OMP END PARALLEL DO

      call halo_exchange(cp,dstruct)
      call halo_exchange(cn,dstruct)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,var)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do var=1,2
          if(aun(var,jedge) > 0._jprw) then
            aun(var,jedge)=aun(var,jedge)*min(limit,cp(var,ip2),cn(var,ip1))
          else
            aun(var,jedge)=aun(var,jedge)*min(limit,cn(var,ip2),cp(var,ip1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_Q


  
  

  subroutine compute_gradient(Q,gradQ,dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(in)    :: Q(:)
    real(kind=jprw), intent(inout) :: gradQ(:,:)
    real(kind=jprw), pointer :: S(:,:)
    real(kind=jprw) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2,jnode
    real(kind=jprw) :: avgQSx(dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)

    ! derivatives 

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(XX,jedge)
      Sy  = S(YY,jedge)
      avgQSx(jedge) = Sx*( Q(ip1) + Q(ip2) )*0.5_jprw
      avgQSy(jedge) = Sy*( Q(ip1) + Q(ip2) )*0.5_jprw
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(XX,jnode) = 0.
      gradQ(YY,jnode) = 0.
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        gradQ(XX,jnode) = gradQ(XX,jnode)+dstruct%sign(jedge,jnode)*avgQSx(iedge)
        gradQ(YY,jnode) = gradQ(YY,jnode)+dstruct%sign(jedge,jnode)*avgQSy(iedge)
      end do
    end do
    !$OMP END PARALLEL DO

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    do jedge = 1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      ip1   = dstruct%edges(iedge,1)
      ip2   = dstruct%edges(iedge,2)
      Sy    = S(YY,iedge)
      avgQ  = ( Q(ip1) + Q(ip2) )*0.5_jprw

      ! correct for wrong Y-derivatives in previous loop,
      gradQ(YY,ip2) = gradQ(YY,ip2) + 2._jprw*Sy*avgQ 
    end do
  end subroutine compute_gradient

  subroutine compute_gradient_abs(Q,gradQ,dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(in)    :: Q(:)
    real(kind=jprw), intent(inout) :: gradQ(:,:)
    real(kind=jprw), pointer :: S(:,:)
    real(kind=jprw) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2,jnode
    real(kind=jprw) :: avgQSx(dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)

    ! derivatives

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(XX,jedge)
      Sy  = S(YY,jedge)
      avgQSx(jedge) = Sx*( abs(Q(ip1)) + abs(Q(ip2)) )*0.5_jprw
      avgQSy(jedge) = Sy*( abs(Q(ip1)) + abs(Q(ip2)) )*0.5_jprw
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(XX,jnode) = 0.
      gradQ(YY,jnode) = 0.
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        gradQ(XX,jnode) = gradQ(XX,jnode)+dstruct%sign(jedge,jnode)*avgQSx(iedge)
        gradQ(YY,jnode) = gradQ(YY,jnode)+dstruct%sign(jedge,jnode)*avgQSy(iedge)
      end do
    end do
    !$OMP END PARALLEL DO

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    do jedge = 1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      ip1   = dstruct%edges(iedge,1)
      ip2   = dstruct%edges(iedge,2)
      Sy    = S(YY,iedge)
      avgQ  = ( Q(ip1) + Q(ip2) )*0.5_jprw

      ! correct for wrong Y-derivatives in previous loop,
      gradQ(YY,ip2) = gradQ(YY,ip2) + 2._jprw*Sy*avgQ
    end do
  end subroutine compute_gradient_abs


  subroutine compute_gradient_tensor(Q,gradQ,dstruct)
    real(kind=jprw), intent(in)  :: Q(:,:)
    real(kind=jprw), intent(out) :: gradQ(:,:)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: S(:,:), D(:)
    real(kind=jprw) :: Sx,Sy
    integer :: jedge,iedge,ip1,ip2,jnode
    real(kind=jprw) :: avgQSx(2,dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(2,dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)
    D   => scalar_field_2d("depth",dstruct)

    ! derivatives

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(XX,jedge)
      Sy  = S(YY,jedge)
      avgQSx(:,jedge) = Sx*( Q(:,ip1) + Q(:,ip2) )*0.5_jprw
      avgQSy(:,jedge) = Sy*( Q(:,ip1) + Q(:,ip2) )*0.5_jprw
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(:,jnode) = 0._jprw
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        gradQ(XXDXX,jnode) = gradQ(XXDXX,jnode)+dstruct%sign(jedge,jnode)*avgQSx(XX,iedge)
        gradQ(XXDYY,jnode) = gradQ(XXDYY,jnode)+dstruct%sign(jedge,jnode)*avgQSy(XX,iedge)
        gradQ(YYDXX,jnode) = gradQ(YYDXX,jnode)+dstruct%sign(jedge,jnode)*avgQSx(YY,iedge)
        gradQ(YYDYY,jnode) = gradQ(YYDYY,jnode)+dstruct%sign(jedge,jnode)*avgQSy(YY,iedge)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine compute_gradient_tensor

  subroutine compute_gradient_tensor_abs(Q,gradQ,dstruct)
    real(kind=jprw), intent(in)  :: Q(:,:)
    real(kind=jprw), intent(out) :: gradQ(:,:)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: S(:,:), D(:)
    real(kind=jprw) :: Sx,Sy
    integer :: jedge,iedge,ip1,ip2,jnode
    real(kind=jprw) :: avgQSx(2,dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(2,dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)
    D   => scalar_field_2d("depth",dstruct)

    ! derivatives

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(XX,jedge)
      Sy  = S(YY,jedge)
      avgQSx(:,jedge) = Sx*( abs(Q(:,ip1)) + abs(Q(:,ip2)) )*0.5_jprw
      avgQSy(:,jedge) = Sy*( abs(Q(:,ip1)) + abs(Q(:,ip2)) )*0.5_jprw
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(:,jnode) = 0._jprw
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        gradQ(XXDXX,jnode) = gradQ(XXDXX,jnode)+dstruct%sign(jedge,jnode)*avgQSx(XX,iedge)
        gradQ(XXDYY,jnode) = gradQ(XXDYY,jnode)+dstruct%sign(jedge,jnode)*avgQSy(XX,iedge)
        gradQ(YYDXX,jnode) = gradQ(YYDXX,jnode)+dstruct%sign(jedge,jnode)*avgQSx(YY,iedge)
        gradQ(YYDYY,jnode) = gradQ(YYDYY,jnode)+dstruct%sign(jedge,jnode)*avgQSy(YY,iedge)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine compute_gradient_tensor_abs

end module mpdata_module
