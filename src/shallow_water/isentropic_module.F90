! =====================================================================
! mpdata3d_module
! This module contains strictly algorithmic subroutines
! - mpdata_gauge       : Non-oscillatory variable sign MPDATA advection
! - compute_gradient   : Compute gradient of a scalar array
! =====================================================================
#include "common/assertions.h"
module mpdata3d_module
  use common_module
  use datastruct_module
  implicit none
  private
  public :: mpdata_D
  public :: mpdata_Q
  public :: compute_gradient
  public :: compute_gradient_tensor

  real(kind=jprw), parameter :: eps = 1.e-12
  integer, parameter, public :: MPDATA_STANDARD=1, MPDATA_GAUGE=2

contains

  subroutine mpdata_D(mpdata_scheme,dt,D,V,VDS,order,limit,dstruct)
    ! For mpdata standard scheme, in the case for positive D,
    !  grad( abs(D) ) == grad( D )
    integer, intent(in) :: mpdata_scheme
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: D(:,:)
    real(kind=jprw), intent(in) :: V(:,:,:)
    real(kind=jprw), intent(in) :: limit
    integer, intent(in) :: order
    real(kind=jprw), intent(out) :: VDS(dstruct%nb_levels,dstruct%nb_edges)

    integer :: jnode, jedge, iedge, jpass, ip1,ip2, jlev
    real(kind=jprw) :: sx, sy, volume_of_two_cells, dDdx, dDdy, Vx, Vy, apos, aneg, Dtmp, D_abs
    real(kind=jprw) :: Dmin(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: Dmax(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: rhin(dstruct%nb_levels)
    real(kind=jprw) :: rhout(dstruct%nb_levels)
    real(kind=jprw) :: cp(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: adv(dstruct%nb_levels)
    real(kind=jprw) :: add
    real(kind=jprw) :: aun(dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: gradD(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: sum_Sabs(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: sum_Dbar(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw), pointer :: vol(:), S(:,:), D0(:,:)
    real(kind=jprw) :: Dbar(2) ! Dabs_bar == Dbar since D > 0 always

    vol     => scalar_field_2d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    D0      => scalar_field_3d("depth_backup",dstruct)

    VDS(:,:) = 0.

    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limit > 0._jprw .and. (order >= 2) ) then
      Dmax(:,:) = -1e10
      Dmin(:,:) =  1e10
      call compute_Dmax_and_Dmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,Sx,Sy,jlev,Vx,Vy,ip1,ip2,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

      Sx = S(XX,jedge)
      Sy = S(YY,jedge)

      do jlev=1,dstruct%nb_levels-1
        Vx = V(XX,jlev,jedge)
        Vy = V(YY,jlev,jedge)

        aun(jlev,jedge) = Vx*Sx + Vy*Sy

        apos = max(0._jprw,aun(jlev,jedge))
        aneg = min(0._jprw,aun(jlev,jedge))
        VDS(jlev,jedge) = D(jlev,ip1)*apos + D(jlev,ip2)*aneg
      end do
    end do
    !$OMP END PARALLEL DO

    !call limit_flux(VDS)

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv(:) = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          add = dstruct%sign(jedge,jnode)
          do jlev=1,dstruct%nb_levels-1
            adv(jlev) = adv(jlev) + add*VDS(jlev,iedge)
          end do
        end do
      end if
      ! Update the unknowns in vertices
      do jlev=1,dstruct%nb_levels-1
        D(jlev,jnode) = max( 0., D(jlev,jnode) - adv(jlev)/vol(jnode) * dt )
      end do
    end do

    !$OMP END PARALLEL DO

    call halo_exchange_3d(D,dstruct) ! Dmax and Dmin could be synced here

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
#if 0
    do jpass=2,order

      ! Compute derivatives for mpdata
      select case (mpdata_scheme)
        case (MPDATA_STANDARD)
          call compute_gradient_abs(D, gradD, dstruct)
        case (MPDATA_GAUGE)
          call compute_gradient(D, gradD, dstruct)
      end select

      call halo_exchange_3d(gradD,dstruct)

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
        call limit_flux(aun)
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

        if( abs(Dtmp) < 1.e-30_jprw ) Dtmp = 0.
        !D(jnode) = Dtmp
        D(jnode) = max(0.,Dtmp)

      enddo
      !$OMP END PARALLEL DO
      call halo_exchange_3d(D,dstruct)


    end do ! other passes
#endif
  contains

    subroutine compute_Dmax_and_Dmin( )
      real(kind=jprw) :: D1, D2
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,D1,D2)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels-1
          D1 = D(jlev,jnode)
          Dmax(jlev,jnode) = max( Dmax(jlev,jnode), D1 )
          Dmin(jlev,jnode) = min( Dmin(jlev,jnode), D1 )
        end do
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          ip2 = dstruct%neighbours(jedge,jnode)
          do jlev=1,dstruct%nb_levels-1
            D2 = D(jlev,ip2)
            Dmax(jlev,jnode) = max( Dmax(jlev,jnode), D2 )
            Dmin(jlev,jnode) = min( Dmin(jlev,jnode), D2 )
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(Dmin,dstruct)
      call halo_exchange_3d(Dmax,dstruct)

    end subroutine compute_Dmax_and_Dmin

    subroutine limit_flux(flux)
      real(kind=jprw), intent(inout) :: flux(:,:)
      real(kind=jprw) :: asignp,asignn

      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,dstruct%nb_nodes
        rhin(:)  = 0.
        rhout(:) = 0.
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          add = dstruct%sign(jedge,jnode)
          do jlev=1,dstruct%nb_levels-1
            apos = max(0._jprw,flux(jlev,iedge))
            aneg = min(0._jprw,flux(jlev,iedge))
            asignp = max(0._jprw,add)
            asignn = min(0._jprw,add)
            rhin(jlev)  = rhin(jlev)  - asignp*aneg - asignn*apos
            rhout(jlev) = rhout(jlev) + asignp*apos + asignn*aneg
          end do
        end do
        do jlev=1,dstruct%nb_levels-1
          cp(jlev,jnode) = ( Dmax(jlev,jnode)-D(jlev,jnode) )*vol(jnode)/( rhin(jlev) * dt + eps )
          cn(jlev,jnode) = ( D(jlev,jnode)-Dmin(jlev,jnode) )*vol(jnode)/( rhout(jlev)* dt + eps )
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(cp,dstruct)
      call halo_exchange_3d(cn,dstruct)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do jlev=1,dstruct%nb_levels-1
          if(flux(jlev,jedge) > 0._jprw) then
            flux(jlev,jedge)=flux(jlev,jedge)*min(limit,cp(jlev,ip2),cn(jlev,ip1))
          else
            flux(jlev,jedge)=flux(jlev,jedge)*min(limit,cn(jlev,ip2),cp(jlev,ip1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO
    end subroutine limit_flux

  end subroutine mpdata_D
  
  
  subroutine mpdata_Q(mpdata_scheme,dt,Q,V,order,limit,dstruct)
    integer, intent(in) :: mpdata_scheme
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: Q(:,:,:)
    real(kind=jprw), intent(in) :: V(:,:,:)
    real(kind=jprw), intent(in) :: limit
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, ip1,ip2, jlev
    real(kind=jprw) :: Sx, Sy, volume_of_two_cells, dQdx(2), dQdy(2), Vx, Vy, Qtmp(2), Q_abs(2), Qbar(2,2)
    real(kind=jprw) :: apos(2), aneg(2), add
    real(kind=jprw) :: Qmin(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: Qmax(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: rhin(2)
    real(kind=jprw) :: rhout(2)
    real(kind=jprw) :: cp(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: cn(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: adv(2,dstruct%nb_levels)
    real(kind=jprw) :: aun(2,dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: fluxv(2,dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: gradQ(4,dstruct%nb_levels,dstruct%nb_nodes)
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
      Qmax(:,:,:) = -1e10
      Qmin(:,:,:) =  1e10
      call compute_Qmax_and_Qmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy,Vx,Vy,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

      Sx = S(XX,jedge)
      Sy = S(YY,jedge)
      do jlev=1,dstruct%nb_levels-1
        Vx = V(XX,jlev,jedge)
        Vy = V(YY,jlev,jedge)

        aun(:,jlev,jedge) = Vx*Sx + Vy*Sy
        apos(:) = max(0._jprw,aun(:,jlev,jedge))
        aneg(:) = min(0._jprw,aun(:,jlev,jedge))
        fluxv(:,jlev,jedge) = Q(:,jlev,ip1)*apos(:) + Q(:,jlev,ip2)*aneg(:)
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv(:,:) = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          add=dstruct%sign(jedge,jnode)
          do jlev=1,dstruct%nb_levels-1
            adv(:,jlev) = adv(:,jlev) + add*fluxv(:,jlev,iedge)
          end do
        end do
      end if
      ! Update the unknowns in vertices
      do jlev=1,dstruct%nb_levels-1
        Q(:,jlev,jnode) = Q(:,jlev,jnode) - adv(:,jlev)/vol(jnode) * dt
      end do
    end do
    !$OMP END PARALLEL DO

    call halo_exchange_3d(Q,dstruct)


    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
#if 0
    do jpass=2,order

      ! Compute derivatives for mpdata
      select case (mpdata_scheme)
        case (MPDATA_STANDARD)
          call compute_gradient_tensor_abs(Q, gradQ, dstruct)
        case (MPDATA_GAUGE)
          call compute_gradient_tensor(Q, gradQ, dstruct)
      end select

      call halo_exchange_3d(gradQ,dstruct)

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
        if (mpdata_scheme == MPDATA_STANDARD) call compute_Qmax_and_Qmin()
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

      call halo_exchange_3d(Q,dstruct)

    end do ! other passes
#endif
  contains

    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1(2), Q2(2)
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,Q1,Q2)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels-1
          Q1(:) = Q(:,jlev,jnode)
          Qmax(:,jlev,jnode) = max( Qmax(:,jlev,jnode), Q1(:) )
          Qmin(:,jlev,jnode) = min( Qmin(:,jlev,jnode), Q1(:) )
        end do
        do jedge = 1,dstruct%nb_neighbours(jnode)
          ip2 = dstruct%neighbours(jedge,jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          do jlev=1,dstruct%nb_levels-1
            Q2(:) = Q(:,jlev,ip2)
            Qmax(:,jlev,jnode) = max( Qmax(:,jlev,jnode), pole_bc(iedge)*Q2(:) )
            Qmin(:,jlev,jnode) = min( Qmin(:,jlev,jnode), pole_bc(iedge)*Q2(:) )
           end do
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(Qmin,dstruct)
      call halo_exchange_3d(Qmax,dstruct)

    end subroutine compute_Qmax_and_Qmin
#if 0
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

      call halo_exchange_3d(cp,dstruct)
      call halo_exchange_3d(cn,dstruct)

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
#endif
  end subroutine mpdata_Q


  

  subroutine compute_gradient(Q,gradQ,dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(in)    :: Q(:,:)
    real(kind=jprw), intent(inout) :: gradQ(:,:,:)
    real(kind=jprw), pointer :: S(:,:), D(:,:)
    real(kind=jprw) :: add
    integer :: jedge,iedge,ip1,ip2,jnode,jlev
    real(kind=jprw) :: avgQS(2,dstruct%nb_levels,dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)
    D   => scalar_field_3d("depth",dstruct)

    ! derivatives 

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      do jlev=1,dstruct%nb_levels
        !avgQS(:,jlev,jedge) = S(:,jedge) * ( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw
        avgQS(XX,jlev,jedge) = S(XX,jedge) * ( D(jlev,ip1)*Q(jlev,ip1) + D(jlev,ip2)*Q(jlev,ip2) ) &
                            & / (  D(jlev,ip1) +  D(jlev,ip2) + eps )
        avgQS(YY,jlev,jedge) = S(YY,jedge) * ( D(jlev,ip1)*Q(jlev,ip1) + D(jlev,ip2)*Q(jlev,ip2) ) &
                            & / (  D(jlev,ip1) +  D(jlev,ip2) + eps )

      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(:,:,jnode) = 0.
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        add = dstruct%sign(jedge,jnode)
        do jlev=1,dstruct%nb_levels
          gradQ(XX,jlev,jnode) = gradQ(XX,jlev,jnode)+add*avgQS(XX,jlev,iedge)
          gradQ(YY,jlev,jnode) = gradQ(YY,jlev,jnode)+add*avgQS(YY,jlev,iedge)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    do jedge = 1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      ip1   = dstruct%edges(iedge,1)
      ip2   = dstruct%edges(iedge,2)
      ! correct for wrong Y-derivatives in previous loop
      do jlev=1,dstruct%nb_levels
        gradQ(YY,jlev,ip2) = gradQ(YY,jlev,ip2) + 2._jprw*avgQS(YY,jlev,iedge)
      end do
    end do

  end subroutine compute_gradient

  subroutine compute_gradient_tensor(Q,gradQ,dstruct)
    real(kind=jprw), intent(in)  :: Q(:,:,:)
    real(kind=jprw), intent(out) :: gradQ(:,:,:)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: S(:,:), D(:,:)
    real(kind=jprw) :: Sx,Sy, add
    integer :: jedge,iedge,ip1,ip2,jnode,jlev
    real(kind=jprw) :: avgQSx(2,dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(2,dstruct%nb_levels,dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)
    D   => scalar_field_3d("depth",dstruct)

    ! derivatives

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(XX,jedge)
      Sy  = S(YY,jedge)
      do jlev=1,dstruct%nb_levels
        !avgQSx(:,jlev,jedge) = Sx*( Q(:,jlev,ip1) + Q(:,jlev,ip2) )*0.5_jprw
        !avgQSy(:,jlev,jedge) = Sy*( Q(:,jlev,ip1) + Q(:,jlev,ip2) )*0.5_jprw

        avgQSx(:,jlev,jedge) = Sx*( D(jlev,ip1)*Q(:,jlev,ip1) + D(jlev,ip2)*Q(:,jlev,ip2) )&
                             & /(D(jlev,ip1)+D(jlev,ip2)+eps)
        avgQSy(:,jlev,jedge) = Sy*( D(jlev,ip1)*Q(:,jlev,ip1) + D(jlev,ip2)*Q(:,jlev,ip2) )&
                             & /(D(jlev,ip1)+D(jlev,ip2)+eps)
      end do

    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(:,:,jnode) = 0._jprw
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        add = dstruct%sign(jedge,jnode)
        do jlev=1,dstruct%nb_levels
          gradQ(XXDXX,jlev,jnode) = gradQ(XXDXX,jlev,jnode)+add*avgQSx(XX,jlev,iedge)
          gradQ(XXDYY,jlev,jnode) = gradQ(XXDYY,jlev,jnode)+add*avgQSy(XX,jlev,iedge)
          gradQ(YYDXX,jlev,jnode) = gradQ(YYDXX,jlev,jnode)+add*avgQSx(YY,jlev,iedge)
          gradQ(YYDYY,jlev,jnode) = gradQ(YYDYY,jlev,jnode)+add*avgQSy(YY,jlev,iedge)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine compute_gradient_tensor


end module mpdata3d_module


! ===================================================================
! isentropic_module
! --------------------
! This module contains subroutines to solve the shallow water eqs
! on a sphere
! - set_time_step : To set the solver's time step
! - setup_isentropic : To modify dual-volumes to map to the sphere
!                         and create necessary fields
! - set_state_zonal_flow : To initialise the state
! - propagate_state : To propagate the state with a given time step
! ===================================================================
module isentropic_module
  use parallel_module
  use common_module
  use datastruct_module
  use mpdata3d_module

  implicit none
  private
  public :: setup_isentropic
  public :: propagate_state
  public :: set_state_zonal_flow
  public :: set_topography
  public :: set_time_step

  real(kind=jprw), parameter :: eps    = 1.e-12
  real(kind=jprw), parameter :: radius = 63.6620e+03
  !real(kind=jprw), parameter :: radius = 6371.22e+03
  real(kind=jprw), parameter :: f0     = 1.4584e-04 !coriolis parameter (=2xearth's omega)
  real(kind=jprw), parameter :: grav   = 9.80616
  real(kind=jprw), parameter :: pi     = acos(-1._jprw)

  real(kind=jprw), parameter :: Rgas   = 287.04 ! gas constant
  real(kind=jprw), parameter :: stf    = 1.02e-05 !stf=bv**2/g
  real(kind=jprw), parameter :: cp     = 3.5*Rgas
  real(kind=jprw), parameter :: cap    = Rgas/cp
  real(kind=jprw), parameter :: pscal  = 1.e5
  real(kind=jprw), parameter :: pscali = 1.e-5
  real(kind=jprw), parameter :: dz     = 100.
  real(kind=jprw), parameter :: pr00   = 1.e5
  real(kind=jprw), parameter :: th00   = 293.15

  real(kind=jprw), parameter :: ibs=1.     ! 0: isentropic/ isosteric --- 1: isopycnic
  real(kind=jprw), parameter :: icp=1.-ibs

  real(kind=jprw), allocatable :: alfs(:), alf0(:), z0(:)

  real(kind=jprw) :: dt_forward = 20.
  integer :: iter = 0

  integer, parameter, public :: EQS_MOMENTUM = 1
  integer, parameter, public :: EQS_VELOCITY = 2

  integer, parameter :: eqs_type = EQS_MOMENTUM

  logical, parameter, public :: STACKED_SHALLOW_WATER = .False.

contains

  ! Helper function
  real(kind=jprw) function alfb(z,z00,stb)
    real(kind=jprw), intent(in) :: z, z00, stb
    alfb = (ibs+icp*th00)*exp(stb*(z-z00))
  end function alfb

  ! Helper function
  real(kind=jprw) function pb(z,z00,stb)
    real(kind=jprw), intent(in) :: z, z00, stb
    pb = ibs*grav/stb*exp(-stb*(z-z00)) + icp*pr00*(1.-grav/(cp*th00*stb) &
       & *(1.-exp(-stb*(z-z00))))**(1./cap)
  end function

  subroutine set_time_step(dt)
    real(kind=jprw), intent(in) :: dt
    dt_forward = dt
  end subroutine set_time_step

  subroutine setup_isentropic(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw) :: y, cos_y, sin_y
    real(kind=jprw), pointer :: coords(:,:), vol(:), hx(:), hy(:), dhxdy_over_G(:), pole_bc(:)
    integer :: jnode, jedge, iedge

    call create_field_in_nodes_2d("hx",1,dstruct)
    call create_field_in_nodes_2d("hy",1,dstruct)
    call create_field_in_nodes_2d("dhxdy_over_G",1,dstruct)
    call create_field_in_edges_2d("pole_bc",1,dstruct)


    coords => vector_field_2d("coordinates",dstruct)
    vol    => scalar_field_2d("dual_volumes",dstruct)
    hx    => scalar_field_2d("hx",dstruct)
    hy    => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)
    !dir$ ivdep
    do jnode=1,dstruct%nb_nodes
      y = coords(YY,jnode)
      cos_y = cos(y)
      sin_y = sin(y)
      hx(jnode) = radius*cos_y
      hy(jnode) = radius
      dhxdy_over_G(jnode) = - sin_y/(radius*max(eps,cos_y))
      vol(jnode) = vol(jnode)*hx(jnode)*hy(jnode)
    enddo

    pole_bc(:) = 1.
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      pole_bc(iedge) = -1.
    end do

    call create_field_in_nodes_2d("coriolis",1,dstruct)
    call create_field_in_nodes_3d("depth",1,dstruct)
    call create_field_in_nodes_3d("velocity",2,dstruct)
    call create_field_in_nodes_3d("forcing",2,dstruct)
    call create_field_in_nodes_3d("depth_backup",1,dstruct)
    call create_field_in_nodes_3d("velocity_backup",2,dstruct)
    call create_field_in_nodes_3d("momentum_backup",2,dstruct)
    call create_field_in_nodes_2d("surface_velocity",2,dstruct)

    call create_field_in_edges_3d("advective_velocity",2,dstruct)
    call create_field_in_nodes_3d("depth_ratio",1,dstruct) ! old/new
    
    call create_field_in_nodes_2d("topography",1,dstruct)
    call create_field_in_nodes_3d("height",1,dstruct)

    call create_field_in_nodes_2d("dlh",1,dstruct)
    call create_field_in_nodes_3d("momentum",2,dstruct)
    call create_field_in_nodes_3d("ambient_depth",1,dstruct)
    call create_field_in_nodes_3d("ambient_momentum",2,dstruct)
    call create_field_in_nodes_3d("pressure",1,dstruct)
    call create_field_in_nodes_3d("montgomery_potential",1,dstruct)
    call create_field_in_nodes_2d("p0",1,dstruct)

    allocate( alfs(dstruct%nb_levels) )
    allocate( alf0(dstruct%nb_levels) )
    allocate( z0(dstruct%nb_levels) )


  end subroutine setup_isentropic

  subroutine propagate_state(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout), target :: dstruct
    real(kind=jprw) :: tend, t0, dt_fwd, tstart
    character(len=200) :: step_info
    real(kind=jprw), pointer :: D(:,:)

    D => scalar_field_3d("depth",dstruct)
    tstart   = dstruct%time
    tend     = dstruct%time+dt
    
    call log_info(" ")
    write(log_str,'(A,I10,A,I10)') "Propagating from time ",int(tstart)," to time ",int(tend); call log_info()
    do while (dstruct%time < tend)
      t0 = dstruct%time
      dt_fwd = min( dt_forward, tend-t0 )
      call step_forward(iter,dt_fwd,dstruct)

      if( log_level <= LOG_LEVEL_INFO ) then
        if( myproc .eq. 0 ) then
          call progress_bar(dstruct%time,tstart,tend)
        end if
      else
        write(step_info,'(A6,I8,A12,F9.1,A12,F8.1,A12,E20.13)') "step = ",iter, &
          & "  time = ",dstruct%time, &
          & "  dt = ",dt_fwd !, "Norm ",L2norm(D)
        CALL LOG_INFO( STEP_INFO )
      end if
    end do

  end subroutine propagate_state




  subroutine step_forward(step,dt,dstruct)
    integer, intent(inout) :: step
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct

    call backup_solution(dstruct)


    if (step == 0) then ! Pre-compute forcing

      call log_debug("Compute forcing")
      call compute_forcing(dstruct)

      call log_debug("Compute advective velocities")
      call compute_advective_velocities(dt,dstruct,"extrapolate")

    end if

    call log_debug("add_forcing_to_solution")
    call add_forcing_to_solution(dt,dstruct)
    
    call log_debug("Advect solution")
    call advect_solution(dt,dstruct)

    call log_debug("Implicit solve")
    if (STACKED_SHALLOW_WATER) then
      call implicit_solve_shallow(dt,dstruct)
    else
      call implicit_solve(dt,dstruct)
    endif

    call log_debug("Compute advective velocities")
    call compute_advective_velocities(dt,dstruct,"advect")

    dstruct%time = dstruct%time+dt
    step = step+1

  end subroutine step_forward



  subroutine set_state_zonal_flow(dstruct)
    type(DataStructure_type), intent(inout)      :: dstruct
    real(kind=jprw), pointer :: D(:,:), cor(:), H0(:), H(:,:), dlh(:), press(:,:), p0(:)
    real(kind=jprw), pointer :: U(:,:,:), coords(:,:), Q(:,:,:), D_amb(:,:), Q_amb(:,:,:)
    real(kind=jprw), pointer :: U0(:,:)
    integer :: jnode, jlev
    real(kind=jprw) :: x,y, zz
    real(kind=jprw), parameter :: v0 = 10.
    real(kind=jprw), parameter :: H00 = grav * 8e3
    real(kind=jprw), parameter :: omega = v0/radius
    real(kind=jprw), parameter :: beta = 0.! pi/4._jprw
    real(kind=jprw), parameter :: pvel = omega*2.

    coords => vector_field_2d("coordinates",dstruct)
    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    U0 => vector_field_2d("surface_velocity",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_3d("height",dstruct)
    dlh => scalar_field_2d("dlh",dstruct)
    D_amb => scalar_field_3d("ambient_depth",dstruct)
    Q_amb => vector_field_3d("ambient_momentum",dstruct)
    press => scalar_field_3d("pressure",dstruct)
    p0 => scalar_field_2d("p0",dstruct)
    Q => vector_field_3d("momentum",dstruct)

    if (STACKED_SHALLOW_WATER) then
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
      do jnode=1,dstruct%nb_nodes
        x=coords(XX,jnode)
        y=coords(YY,jnode)
        cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
        H(:,jnode)     = (H00-radius**2*(f0+pvel)*0.5*pvel &
                     & *(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2)/grav
        D(:,jnode)     = max(0., H(:,jnode) - H0(jnode) )
        U(XX,:,jnode)  =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y)
        U(YY,:,jnode)  = -pvel*sin(x)*sin(beta)*radius
        H(:,jnode) = H0(jnode) + D(:,jnode)
        do jlev=1,dstruct%nb_levels
          Q(:,jlev,jnode) = D(jlev,jnode) * U(:,jlev,jnode)
        end do
      end do
      !$OMP END PARALLEL DO
      call compute_forcing_shallow(dstruct)

    else
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
      do jnode=1,dstruct%nb_nodes
        x=coords(XX,jnode)
        y=coords(YY,jnode)
        cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
        U0(XX,jnode) =  omega*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*max(0._jprw,cos(y))
        U0(YY,jnode) = -omega*sin(x)*sin(beta)*radius

        dlh(jnode)   = -radius**2*(f0+omega)*0.5_jprw*omega/grav &
                     & *(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2

        ! Initial Z at levels          This is where pressure and height are defined
        do jlev=1,dstruct%nb_levels
          Z0(jlev) = (jlev-1)*dz
          alf0(jlev) = alfb(Z0(jlev),0._jprw,stf)
        end do
        ! Z at half levels (staggered)  This is where unknowns D,Q are defined
        do jlev=1,dstruct%nb_levels-1
          alfs(jlev) = alfb(Z0(jlev)+dz*0.5_jprw,0._jprw,stf)
        end do
        ! TODO call absorbers(taui,z0,nlm)

        ! Set height
        do jlev=1,dstruct%nb_levels
          H(jlev,jnode) = max(H0(jnode), Z0(jlev)+dlh(jnode))
        end do

        ! At top boundary
        zz = max( Z0(dstruct%nb_levels), H0(jnode) )
        press(dstruct%nb_levels,jnode) = pb(zz,0._jprw,stf)
        p0(jnode) = press(dstruct%nb_levels,jnode)

        ! Convert pressure to exner in case ibs==0 (isentropic)
        press(dstruct%nb_levels,jnode) = ibs*p0(jnode)+icp*cp*(p0(jnode)/pr00)**cap

        ! Integrate hydrostatic relation
        do jlev=dstruct%nb_levels-1,1,-1
           press(jlev,jnode) = press(jlev+1,jnode)+grav*(H(jlev+1,jnode)-H(jlev,jnode))/alfs(jlev)
        enddo

        ! Convert exner to pressure in case ibs==0
        do jlev=1,dstruct%nb_levels
          press(jlev,jnode) = ibs*press(jlev,jnode)+icp*pr00*(press(jlev,jnode)/cp)**(1/cap)
        enddo

        ! Set ambient and initial state
        do jlev=1,dstruct%nb_levels-1
          D_amb(jlev,jnode)    = (press(jlev,jnode)-press(jlev+1,jnode))*pscali
          Q_amb(XX,jlev,jnode) = (U0(XX,jnode)-0.01*cos(y)*H(jlev,jnode))*D_amb(jlev,jnode)
          Q_amb(YY,jlev,jnode) = U0(YY,jnode)*D_amb(jlev,jnode)
          D(jlev,jnode)   = D_amb(jlev,jnode)
          Q(:,jlev,jnode) = Q_amb(:,jlev,jnode)
          U(:,jlev,jnode) = Q(:,jlev,jnode)/max(eps,D(jlev,jnode))
        enddo
      end do
      !$OMP END PARALLEL DO
      call compute_forcing(dstruct)

    end if

  end subroutine set_state_zonal_flow


  subroutine set_topography(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:), pointer :: H0
    real(kind=jprw), dimension(:,:), pointer :: coords
    real(kind=jprw) :: zlatc, zlonc, zr, zz, rad, x, y, amp
    integer :: jnode

    if (STACKED_SHALLOW_WATER) then
      call set_topography_mountain(800._jprw,dstruct)
    else
      zlatc=0.
      zlonc=1.5*pi
      zr = 5000.
      amp = 50.

      H0 => scalar_field_2d("topography",dstruct)
      coords => vector_field_2d("coordinates",dstruct)

      do jnode=1,dstruct%nb_nodes
        x = coords(XX,jnode)
        y = coords(YY,jnode)
        zz = sin(zlatc)*sin(y)+cos(zlatc)*cos(y)*cos(x-zlonc)
        rad = radius*acos(zz) / zr
        H0(jnode) = amp/(sqrt(1.+rad**2))**3
      end do
    end if
  end subroutine set_topography


  subroutine set_topography_mountain(amplitude,dstruct)
    real(kind=jprw), intent(in) :: amplitude ! amplitude of hill
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:), pointer :: H0
    real(kind=jprw), dimension(:,:), pointer :: coords
    real(kind=jprw) :: rad = 2.*pi/18. ! radius of hill
    real(kind=jprw) :: xcent = 3.*pi/2.  ! centre of hill
    real(kind=jprw) :: ycent = pi/6.*1.
    real(kind=jprw) :: gamm = 0.5 ! slope of hill
    real(kind=jprw) :: dist, xlon, ylat
    integer :: jnode

    H0 => scalar_field_2d("topography",dstruct)
    coords => vector_field_2d("coordinates",dstruct)

    do jnode=1,dstruct%nb_nodes
      xlon = coords(XX,jnode)
      ylat = coords(YY,jnode)

      dist = 2.*sqrt( (cos(ylat)*sin( (xlon-xcent)/2 ) )**2 &
        &     + sin((ylat-ycent)/2)**2 )

      H0(jnode) = max(0.,amplitude * (1.-gamm*dist/rad))
    end do
  end subroutine set_topography_mountain


  subroutine backup_solution(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: D(:,:), D0(:,:)
    real(kind=jprw), pointer :: Q(:,:,:), Q0(:,:,:), U(:,:,:), U0(:,:,:)
    integer :: jnode
    
    D  => scalar_field_3d("depth",dstruct)
    D0 => scalar_field_3d("depth_backup",dstruct)
    Q  => vector_field_3d("momentum",dstruct)
    Q0 => vector_field_3d("momentum_backup",dstruct)
    U  => vector_field_3d("velocity",dstruct)
    U0 => vector_field_3d("velocity_backup",dstruct)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      D0(:,jnode) = D(:,jnode)
      Q0(:,:,jnode) = Q(:,:,jnode)
      U0(:,:,jnode) = U(:,:,jnode)
    end do
    !$OMP END PARALLEL DO
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
    integer :: jnode, jedge, iedge, jlev, ip1, ip2
    real(kind=jprw), dimension(:),     pointer :: hx, hy, vol
    real(kind=jprw), dimension(:,:),   pointer :: D, D0, coords
    real(kind=jprw), dimension(:,:,:), pointer :: U, U0, R, Vedges, Q, Q0
    real(kind=jprw) :: Vnodes(2,dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: grad_Vnodes(4,dstruct%nb_levels,dstruct%nb_nodes)

    coords => vector_field_2d("coordinates",dstruct)
    Vedges => vector_field_3d("advective_velocity",dstruct)
    D      => scalar_field_3d("depth",dstruct)
    D0     => scalar_field_3d("depth_backup",dstruct)
    Q      => vector_field_3d("momentum",dstruct)
    Q0     => vector_field_3d("momentum_backup",dstruct)
    U      => vector_field_3d("velocity",dstruct)
    U0     => vector_field_3d("velocity_backup",dstruct)
    hx     => scalar_field_2d("hx",dstruct)
    hy     => scalar_field_2d("hy",dstruct)
    R      => vector_field_3d("forcing",dstruct)
    vol    => scalar_field_2d("dual_volumes",dstruct)

    if( option .eq. "advect") then
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Dpos)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels
          Dpos = max(eps, D(jlev,jnode))
          Vnodes(:,jlev,jnode) = Q(:,jlev,jnode)/Dpos
          U(:,jlev,jnode)      = Q(:,jlev,jnode)/Dpos
        end do
      end do
      !$OMP END PARALLEL DO

      call compute_gradient_tensor( Vnodes, grad_Vnodes, dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Vx,Vy,dVxdx,dVxdy,dVydx,dVydy)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels
          Ux = U(XX,jlev,jnode)
          Uy = U(YY,jlev,jnode)
          Vx = Ux
          Vy = Uy
          dVxdx = grad_Vnodes(XXDXX,jlev,jnode)*hy(jnode)/vol(jnode)
          dVxdy = grad_Vnodes(XXDYY,jlev,jnode)*hx(jnode)/vol(jnode)
          dVydx = grad_Vnodes(YYDXX,jlev,jnode)*hy(jnode)/vol(jnode)
          dVydy = grad_Vnodes(YYDYY,jlev,jnode)*hx(jnode)/vol(jnode)
          Dpos = max(eps,D(jlev,jnode))
          Vnodes(XX,jlev,jnode) = ( Ux - 0.5*dt*(Vx*dVxdx+Vy*dVxdy) + 0.5*dt*R(XX,jlev,jnode)/Dpos ) * hy(jnode)
          Vnodes(YY,jlev,jnode) = ( Uy - 0.5*dt*(Vx*dVydx+Vy*dVydy) + 0.5*dt*R(YY,jlev,jnode)/Dpos ) * hx(jnode)
          !Vnodes(XX,jlev,jnode) = 10.*cos(coords(YY,jnode)) * hy(jnode)
          !Vnodes(YY,jlev,jnode) = 0.  * hx(jnode)
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d( Vnodes, dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,jlev)
      do jedge=1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do jlev=1,dstruct%nb_levels
          !Vedges(:,jlev,jedge) = (Vnodes(:,jlev,ip1)+Vnodes(:,jlev,ip2))*0.5_jprw
          Vedges(:,jlev,jedge) = (D(jlev,ip1)*Vnodes(:,jlev,ip1)+D(jlev,ip2)*Vnodes(:,jlev,ip2)) &
                               & /(D(jlev,ip1)+D(jlev,ip2)+eps)
        end do
      end do
      !$OMP END PARALLEL DO

    else if( option .eq. "extrapolate") then
      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Dpos,D0pos,Ux,Uy,U0x,U0y)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels-1
          Dpos  = max(eps, D(jlev,jnode))
          D0pos  = max(eps, D0(jlev,jnode))
          Ux    = Q(XX,jlev,jnode)/Dpos
          Uy    = Q(YY,jlev,jnode)/Dpos
          U0x   = Q0(XX,jlev,jnode)/D0pos
          U0y   = Q0(YY,jlev,jnode)/D0pos
          Vnodes(XX,jlev,jnode) = ( 1.5_jprw*Ux - 0.5_jprw*U0x ) * hy(jnode)
          Vnodes(YY,jlev,jnode) = ( 1.5_jprw*Uy - 0.5_jprw*U0y ) * hx(jnode)
        end do
      end do
      !$OMP END PARALLEL DO

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge=1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do jlev=1,dstruct%nb_levels-1
          !Vedges(:,jlev,jedge) = (Vnodes(:,jlev,ip1)+Vnodes(:,jlev,ip2))*0.5_jprw
          Vedges(:,jlev,jedge) = (D(jlev,ip1)*Vnodes(:,jlev,ip1)+D(jlev,ip2)*Vnodes(:,jlev,ip2)) &
                               & /(D(jlev,ip1)+D(jlev,ip2)+eps)
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      do jlev=1,dstruct%nb_levels-1
        Vedges(YY,jlev,iedge) = 0.
      end do
    end do

  end subroutine compute_advective_velocities

  subroutine compute_forcing_shallow(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode,jlev
    real(kind=jprw) :: Ux, Uy, Qx, Qy, Dpos
    real(kind=jprw), dimension(:),   pointer :: H0, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:),   pointer :: H, D
    real(kind=jprw), dimension(:,:,:), pointer :: U, Q, R
    real(kind=jprw) :: grad_H(2,dstruct%nb_levels,dstruct%nb_nodes)
    vol => scalar_field_2d("dual_volumes",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    H => scalar_field_3d("height",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    Q => vector_field_3d("momentum",dstruct)
    R => vector_field_3d("forcing",dstruct)

    do jnode=1,dstruct%nb_nodes
      H(:,jnode) = H0(jnode) + D(:,jnode)
    end do
    call compute_gradient( H, grad_H, dstruct )
    call halo_exchange_3d(grad_H,dstruct)

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
    do jnode=1,dstruct%nb_nodes
      do jlev=1,dstruct%nb_levels
        Qx    = Q(XX,jlev,jnode)
        Qy    = Q(YY,jlev,jnode)
        Dpos  = max(eps, D(jlev,jnode))
        R(XX,jlev,jnode) = -grav*D(jlev,jnode)*grad_H(XX,jlev,jnode)*hy(jnode)/vol(jnode) &
          &           + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/Dpos
        R(YY,jlev,jnode) = -grav*D(jlev,jnode)*grad_H(YY,jlev,jnode)*hx(jnode)/vol(jnode) &
          &           - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/Dpos
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine compute_forcing_shallow

  subroutine compute_forcing(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, jlev
    real(kind=jprw) :: Qx, Qy
    real(kind=jprw), dimension(:),   pointer :: H0, hx, hy, cor, p0, vol, dhxdy_over_G
    real(kind=jprw), dimension(:,:), pointer :: D, coords, M, press, H
    real(kind=jprw), dimension(:,:,:), pointer :: Q, R
    real(kind=jprw) :: grad_M(2,dstruct%nb_levels,dstruct%nb_nodes)

    coords => vector_field_2d("coordinates",dstruct)
    vol => scalar_field_2d("dual_volumes",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    H0 => scalar_field_2d("topography",dstruct)
    D  => scalar_field_3d("depth",dstruct)
    H  => scalar_field_3d("height",dstruct)
    Q  => vector_field_3d("momentum",dstruct)
    R  => vector_field_3d("forcing",dstruct)
    M  => scalar_field_3d("montgomery_potential",dstruct)
    press  => scalar_field_3d("pressure",dstruct)
    p0 => scalar_field_2d("p0",dstruct)

    do jnode=1,dstruct%nb_nodes

      if( STACKED_SHALLOW_WATER ) then
        do jlev=1,dstruct%nb_levels
          H(jlev,jnode) = H0(jnode) + D(jlev,jnode)
          M(jlev,jnode) = H(jlev,jnode)*grav
        end do
      else
        ! Get pressure from pressure thickness
        press(dstruct%nb_levels,jnode) = p0(jnode)
        do jlev=dstruct%nb_levels-1,1,-1
          press(jlev,jnode)= press(jlev+1,jnode)+D(jlev,jnode)*pscal
        enddo

        ! Convert pressure to exner in case ibs==0 (isentropic)
        do jlev=1,dstruct%nb_levels
          press(jlev,jnode) = ibs*press(jlev,jnode)+icp*cp*(press(jlev,jnode)/pr00)**cap
        end do

        ! Define Montegomery potential (alias Bernoulli head)
        M(1,jnode) = alf0(1)*press(1,jnode)+grav*H0(jnode) + (alfs(1)-alf0(1))*press(1,jnode)

        do jlev=2,dstruct%nb_levels-1
           M(jlev,jnode) = M(jlev-1,jnode) + (alfs(jlev)-alf0(jlev-1))*press(jlev,jnode)
        end do

        ! Convert exner to pressure in case ibs==0
        do jlev=1,dstruct%nb_levels
          press(jlev,jnode) = ibs*press(jlev,jnode)+icp*pr00*(press(jlev,jnode)/cp)**(1/cap)
        end do
      end if
    end do


    call compute_gradient( M, grad_M, dstruct )
    call halo_exchange_3d( grad_M, dstruct)

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
    do jnode=1,dstruct%nb_nodes
      do jlev=1,dstruct%nb_levels
        Qx    = Q(XX,jlev,jnode)
        Qy    = Q(YY,jlev,jnode)
        R(XX,jlev,jnode) = -D(jlev,jnode)*grad_M(XX,jlev,jnode)*hy(jnode)/vol(jnode) &
          &              + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/max(eps,D(jlev,jnode))
        R(YY,jlev,jnode) = -D(jlev,jnode)*grad_M(YY,jlev,jnode)*hx(jnode)/vol(jnode) &
          &              - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/max(eps,D(jlev,jnode))
      end do
    end do
    !$OMP END PARALLEL DO

    !call halo_exchange_3d(R,dstruct)

  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: U(:,:,:) , Q(:,:,:), R(:,:,:), D(:,:)
    integer :: jnode

    D => scalar_field_3d("depth",dstruct)
    Q => vector_field_3d("momentum",dstruct)
    U => vector_field_3d("velocity",dstruct)
    R => vector_field_3d("forcing",dstruct)
    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      Q(:,:,jnode) = Q(:,:,jnode) + 0.5_jprw*dt*R(:,:,jnode)
    end do
    !$OMP END PARALLEL DO
  end subroutine add_forcing_to_solution


  subroutine implicit_solve_shallow(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, m, jlev
    real(kind=jprw) :: Ux, Uy, Qx, Qy, Rx, Ry, Dpos
    real(kind=jprw) :: Ux_adv, Uy_adv, Qx_adv, Qy_adv, Rx_exp, Ry_exp

    real(kind=jprw), dimension(:),   pointer   :: H0, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:),   pointer :: H, D
    real(kind=jprw), dimension(:,:,:), pointer :: U, Q, R
    real(kind=jprw) :: grad_H(2,dstruct%nb_levels,dstruct%nb_nodes)

    vol => scalar_field_2d("dual_volumes",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_3d("height",dstruct)
    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    Q => vector_field_3d("momentum",dstruct)
    R => vector_field_3d("forcing",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    do jnode=1,dstruct%nb_nodes
      H(:,jnode) = H0(jnode) + D(:,jnode)
    end do
    call compute_gradient( H, grad_H, dstruct )
    call halo_exchange_3d( grad_H, dstruct )

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Qx,Qy,Dpos,Rx_exp,Ry_exp,Qx_adv,Qy_adv,m,Rx,Ry)
    do jnode=1,dstruct%nb_nodes
      do jlev=1,dstruct%nb_levels
        Dpos  = D(jlev,jnode)

        Qx    = Q(XX,jlev,jnode)
        Qy    = Q(YY,jlev,jnode)

        Qx_adv = Qx
        Qy_adv = Qy

        Rx_exp = -grav*D(jlev,jnode)*grad_H(XX,jlev,jnode)*hy(jnode)/vol(jnode)
        Ry_exp = -grav*D(jlev,jnode)*grad_H(YY,jlev,jnode)*hx(jnode)/vol(jnode)

        do m=1,3 ! Three iterations at most is enough to converge
          Rx = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/D(jlev,jnode)
          Ry = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/D(jlev,jnode)
          Qx = Qx_adv + 0.5_jprw*dt*Rx
          Qy = Qy_adv + 0.5_jprw*dt*Ry
        end do
        Q(XX,jlev,jnode) = Qx
        Q(YY,jlev,jnode) = Qy
        R(XX,jlev,jnode) = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/D(jlev,jnode)
        R(YY,jlev,jnode) = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/D(jlev,jnode)
        U(:,jlev,jnode)  = Q(:,jlev,jnode) / D(jlev,jnode)
      end do
    end do
    !$OMP END PARALLEL DO

  end subroutine implicit_solve_shallow


  subroutine implicit_solve(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, jlev, m
    real(kind=jprw), dimension(:),   pointer :: H0, hx, hy, cor, p0, vol, dhxdy_over_G
    real(kind=jprw), dimension(:,:), pointer :: D, coords, Mont, press, H
    real(kind=jprw), dimension(:,:,:), pointer :: Q, R
    real(kind=jprw) :: grad_M(2,dstruct%nb_levels,dstruct%nb_nodes)

    real(kind=jprw) :: Qx, Qy, Rx, Ry
    real(kind=jprw) :: Qx_adv, Qy_adv, Rx_exp, Ry_exp

    coords => vector_field_2d("coordinates",dstruct)
    vol => scalar_field_2d("dual_volumes",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    H0 => scalar_field_2d("topography",dstruct)
    D  => scalar_field_3d("depth",dstruct)
    H  => scalar_field_3d("height",dstruct)
    Q  => vector_field_3d("momentum",dstruct)
    R  => vector_field_3d("forcing",dstruct)
    Mont  => scalar_field_3d("montgomery_potential",dstruct)
    press  => scalar_field_3d("pressure",dstruct)
    p0 => scalar_field_2d("p0",dstruct)

    do jnode=1,dstruct%nb_nodes
    if( STACKED_SHALLOW_WATER ) then
      do jlev=1,dstruct%nb_levels-1
        H(jlev,jnode) = H0(jnode) + D(jlev,jnode)
        Mont(jlev,jnode) = H(jlev,jnode)*grav
      end do
    else

      ! Get pressure from pressure thickness
      press(dstruct%nb_levels,jnode) = p0(jnode)
      do jlev=dstruct%nb_levels-1,1,-1
        press(jlev,jnode) = press(jlev+1,jnode)+D(jlev,jnode)*pscal
      enddo

      ! Convert pressure to exner in case ibs==0 (isentropic)
      do jlev=1,dstruct%nb_levels
        press(jlev,jnode) = ibs*press(jlev,jnode)+icp*cp*(press(jlev,jnode)/pr00)**cap
      end do

      ! Define Montegomery potential (alias Bernoulli head)
      Mont(1,jnode) = alf0(1)*press(1,jnode)+grav*H0(jnode) + (alfs(1)-alf0(1))*press(1,jnode)

      do jlev=2,dstruct%nb_levels-1
         Mont(jlev,jnode) = Mont(jlev-1,jnode) + (alfs(jlev)-alf0(jlev-1))*press(jlev,jnode)
      end do

      ! Convert exner to pressure in case ibs==0
      do jlev=1,dstruct%nb_levels
        press(jlev,jnode) = ibs*press(jlev,jnode)+icp*pr00*(press(jlev,jnode)/cp)**(1/cap)
      end do

      H(1,jnode) = H0(jnode)
      do jlev=2,dstruct%nb_levels
        H(jlev,jnode)=H(jlev-1,jnode) -alfs(jlev-1)*(press(jlev,jnode)-press(jlev-1,jnode))/grav
      enddo

      end if
    end do

    call compute_gradient( Mont, grad_M, dstruct )
    call halo_exchange_3d(grad_M,dstruct)

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
    do jnode=1,dstruct%nb_nodes
      do jlev=1,dstruct%nb_levels-1
        Qx    = Q(XX,jlev,jnode)
        Qy    = Q(YY,jlev,jnode)
        Rx_exp = -D(jlev,jnode)*grad_M(XX,jlev,jnode)*hy(jnode)/vol(jnode)
        Ry_exp = -D(jlev,jnode)*grad_M(YY,jlev,jnode)*hx(jnode)/vol(jnode)
        Qx_adv = Qx
        Qy_adv = Qy
        do m=1,3 ! Three iterations at most is enough to converge
          Rx = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/max(eps,D(jlev,jnode))
          Ry = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/max(eps,D(jlev,jnode))
          Qx = Qx_adv + 0.5_jprw*dt*Rx
          Qy = Qy_adv + 0.5_jprw*dt*Ry
        end do
        Q(XX,jlev,jnode) = Qx
        Q(YY,jlev,jnode) = Qy
        R(XX,jlev,jnode) = Rx_exp + cor(jnode)*Qy - dhxdy_over_G(jnode)*Qx*Qy/max(eps,D(jlev,jnode))
        R(YY,jlev,jnode) = Ry_exp - cor(jnode)*Qx + dhxdy_over_G(jnode)*Qx*Qx/max(eps,D(jlev,jnode))
      end do
    end do
    !$OMP END PARALLEL DO

!    call halo_exchange_3d(Q,dstruct)
!    call halo_exchange_3d(R,dstruct)

  end subroutine implicit_solve



  subroutine advect_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:,:),   pointer :: D, D0, DR
    real(kind=jprw), dimension(:,:,:), pointer :: Q, V
    real(kind=jprw) :: VDS(dstruct%nb_levels,dstruct%nb_edges)
    integer :: jnode
    
    D => scalar_field_3d("depth",dstruct)
    D0 => scalar_field_3d("depth_backup",dstruct)
    Q => vector_field_3d("momentum",dstruct)
    V => vector_field_3d("advective_velocity",dstruct)
   
    !    mpdata_D( scheme,       time, variable, velocity, VDS,  order, limit,   dstruct )
    call mpdata_D( MPDATA_GAUGE, dt,   D,        V,        VDS,  1,     1._jprw, dstruct )
    call log_debug("advected D")
    !    mpdata_Q( scheme,       time, variable, velocity,  order, limit,   dstruct )
    call mpdata_Q( MPDATA_GAUGE, dt,   Q,        V,         1,     1._jprw, dstruct )
    call log_debug("advected Q")

  end subroutine advect_solution

end module isentropic_module
