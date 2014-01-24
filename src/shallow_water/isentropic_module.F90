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
  public :: mpdata_gauge_D
  public :: mpdata_gauge_Q
  public :: compute_gradient

  real(kind=jprw), parameter :: eps = 1.e-6

contains

  subroutine mpdata_gauge_D(dt,Q,V,VDS,order,limited,Q_is_vector_component,dstruct )
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: Q(:,:)
    real(kind=jprw), intent(in) :: V(:,:,:)
    logical, intent(in) :: Q_is_vector_component, limited
    integer, intent(in) :: order
    real(kind=jprw), intent(out) :: VDS(dstruct%nb_levels,dstruct%nb_edges)
    
    integer :: jnode, jedge, iedge, jpass, ip1,ip2, jlev
    real(kind=jprw) :: sx, sy, volume_of_two_cells, dQdx, dQdy, Vx, Vy, apos, aneg
    real(kind=jprw) :: Qmin(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: Qmax(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: rhin(dstruct%nb_levels)
    real(kind=jprw) :: rhout(dstruct%nb_levels)
    real(kind=jprw) :: cp(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: fluxv(dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: gradQ(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw), pointer :: vol(:,:), S(:,:), pole_bc(:)
    

    vol     => scalar_field_3d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)

    VDS(:,:) = 0.

    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limited .and. (order .ge. 2) ) then
      Qmax(:,:) = -1e10
      Qmin(:,:) =  1e10
      call log_debug("Compute Qmax and Qmin")
      call compute_Qmax_and_Qmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,jlev,Sx,Sy,Vx,Vy,ip1,ip2,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      Sx = S(jedge,XX)
      Sy = S(jedge,YY)
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

      do jlev=1,dstruct%nb_levels
        Vx = V(jlev,jedge,XX)
        Vy = V(jlev,jedge,YY)

        aun(jlev,jedge) = Vx*Sx + Vy*Sy 

        apos = max(0._jprw,aun(jlev,jedge))
        aneg = min(0._jprw,aun(jlev,jedge))
        fluxv(jlev,jedge) = Q(jlev,ip1)*apos + Q(jlev,ip2)*aneg
        VDS(jlev,jedge) = fluxv(jlev,jedge)
      end do
    enddo
    !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge,jlev)
    do jnode=1,dstruct%nb_nodes
      adv = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          do jlev=1,dstruct%nb_levels
            adv = adv + dstruct%sign(jedge,jnode)*fluxv(jlev,iedge)
          end do
        enddo
      endif
      ! Update the unknowns in vertices
      do jlev=1,dstruct%nb_levels
        Q(jlev,jnode) = Q(jlev,jnode) - adv/vol(jlev,jnode) * dt
      end do
    enddo
    !$OMP END PARALLEL DO

    call halo_exchange_3d(Q,dstruct) ! Qmax and Qmin could be synced here

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector_component, dstruct)

      call halo_exchange_3d(gradQ,dstruct)

      ! Compute antidiffusive normal velocity in faces

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,jlev,volume_of_two_cells,dQdx,dQdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        do jlev=1,dstruct%nb_levels
          ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
          volume_of_two_cells = vol(jlev,ip1) + vol(jlev,ip2)
          dQdx = (gradQ(jlev,ip1,XX)+gradQ(jlev,ip2,XX)) / volume_of_two_cells
          dQdy = (gradQ(jlev,ip1,YY)+gradQ(jlev,ip2,YY)) / volume_of_two_cells
          Vx = V(jlev,jedge,XX)
          Vy = V(jlev,jedge,YY)
          ! variable sign option with asymptotic analysis, (mpdata gauge)
          aun(jlev,jedge) = abs(aun(jlev,jedge))*(Q(jlev,ip2)-Q(jlev,ip1))*0.5_jprw &
            &          -0.5_jprw*dt*aun(jlev,jedge)*(Vx*dQdx+Vy*dQdy)
        end do
      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limited) then
        call limit_antidiffusive_velocity()
      endif

      do jedge=1,dstruct%nb_edges
        do jlev=1,dstruct%nb_levels
          VDS(jlev,jedge) = VDS(jlev,jedge) + aun(jlev,jedge)
        end do
      end do

      ! Compute fluxes from (limited) antidiffusive velocity
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge,jlev)
      do jnode=1,dstruct%nb_nodes
        adv = 0.0
        if(dstruct%nb_neighbours(jnode) > 1) then
          do jedge = 1,dstruct%nb_neighbours(jnode)
            iedge = dstruct%my_edges(jedge,jnode)
            do jlev=1,dstruct%nb_levels
              fluxv(jlev,iedge) = aun(jlev,iedge)
              adv = adv + dstruct%sign(jedge,jnode)*fluxv(jlev,iedge)
            end do
          enddo
        endif
        ! Update the unknowns in vertices
        do jlev=1,dstruct%nb_levels
          Q(jlev,jnode) = max( Q(jlev,jnode) - adv/vol(jlev,jnode) * dt, 0. )
        end do
      enddo
      !$OMP END PARALLEL DO
      call halo_exchange_3d(Q,dstruct)

    end do ! other passes

  contains
    
    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1, Q2    
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,Q1,Q2)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels
          Q1 = Q(jlev,jnode)
          Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q1 )
          Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q1 )
        end do
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          do jlev=1,dstruct%nb_levels
            Q2 = Q(jlev,dstruct%neighbours(jedge,jnode))
            if (.not. Q_is_vector_component) then
              Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q2 )
              Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q2 )
            else
              Qmax(jlev,jnode) = max( Qmax(jlev,jnode), pole_bc(iedge)*Q2 )
              Qmin(jlev,jnode) = min( Qmin(jlev,jnode), pole_bc(iedge)*Q2 )
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(Qmin,dstruct)
      call halo_exchange_3d(Qmax,dstruct)

    end subroutine compute_Qmax_and_Qmin

    subroutine limit_antidiffusive_velocity

      real(kind=jprw) :: asignp,asignn
      real(kind=jprw) :: limit = 1.  ! 1: second order, 0: first order

      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,dstruct%nb_nodes
        rhin  = 0.
        rhout = 0.
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          do jlev=1,dstruct%nb_levels
            apos = max(0._jprw,aun(jlev,iedge))
            aneg = min(0._jprw,aun(jlev,iedge))
            rhin(jlev)  = rhin(jlev)  - asignp*aneg - asignn*apos
            rhout(jlev) = rhout(jlev) + asignp*apos + asignn*aneg
          end do
        end do
        do jlev=1,dstruct%nb_levels
          cp(jlev,jnode) = ( Qmax(jlev,jnode)-Q(jlev,jnode) )*vol(jlev,jnode)/( rhin(jlev) * dt + eps )
          cn(jlev,jnode) = ( Q(jlev,jnode)-Qmin(jlev,jnode) )*vol(jlev,jnode)/( rhout(jlev)* dt + eps )
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(cp,dstruct)
      call halo_exchange_3d(cn,dstruct)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do jlev=1,dstruct%nb_levels
          if(aun(jlev,jedge) > 0._jprw) then
            aun(jlev,jedge)=aun(jlev,jedge)*min(limit,cp(jlev,ip2),cn(jlev,ip1))
          else
            aun(jlev,jedge)=aun(jlev,jedge)*min(limit,cn(jlev,ip2),cp(jlev,ip1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_gauge_D
  
  
  
  
  
  
  
  subroutine mpdata_gauge_Q(dt,Q,VDS,DR,D,order,limited,Q_is_vector_component,dstruct)
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: Q(:,:), D(:,:)
    real(kind=jprw), intent(in) :: VDS(:,:), DR(:,:)
    logical, intent(in) :: Q_is_vector_component, limited
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, jlev, ip1,ip2
    real(kind=jprw) :: Sx, Sy, Ssqr, volume_of_two_cells, dQdx, dQdy, Vx, Vy
    real(kind=jprw) :: apos, aneg, x1, x2, y1, y2, length
    real(kind=jprw) :: Qmin(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: Qmax(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: rhin(dstruct%nb_levels)
    real(kind=jprw) :: rhout(dstruct%nb_levels)
    real(kind=jprw) :: cp(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: fluxv(dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: gradQ(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw) :: V(dstruct%nb_levels,dstruct%nb_edges,2)
    real(kind=jprw) :: VDnodes(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw) :: Lengths(dstruct%nb_nodes)
    real(kind=jprw) :: Qtmp(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: volD(dstruct%nb_levels,dstruct%nb_nodes)

    real(kind=jprw), pointer :: vol(:,:), S(:,:), pole_bc(:), coords(:,:)

    vol     => scalar_field_3d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)
    coords  => vector_field_2d("coordinates",dstruct)

    volD(:,:) = vol(:,:)*D(:,:)
    Lengths(:) = 0.
    VDnodes(:,:,:) = 0.
    do jedge=1,dstruct%nb_edges
      Sx = S(jedge,XX)
      Sy = S(jedge,YY)
      Ssqr =  Sx*Sx + Sy*Sy 
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      x1 = coords(ip1,XX)
      x2 = coords(ip2,XX)
      y1 = coords(ip1,YY)
      y2 = coords(ip2,YY)
      length = 1._jprw/sqrt( (x2-x1)**2 + (y2-y1)**2 )
      Lengths(ip1) = Lengths(ip1) + length
      Lengths(ip2) = Lengths(ip2) + length
  
      do jlev=1,dstruct%nb_levels
        V(jlev,jedge,XX) = Sx/Ssqr * VDS(jlev,jedge)
        V(jlev,jedge,YY) = Sy/Ssqr * VDS(jlev,jedge)

        VDnodes(jlev,ip1,XX) = VDnodes(jlev,ip1,XX) + V(jlev,jedge,XX) * length
        VDnodes(jlev,ip2,XX) = VDnodes(jlev,ip2,XX) + V(jlev,jedge,XX) * length
        VDnodes(jlev,ip1,YY) = VDnodes(jlev,ip1,YY) + V(jlev,jedge,YY) * length
        VDnodes(jlev,ip2,YY) = VDnodes(jlev,ip2,YY) + V(jlev,jedge,YY) * length
      end do
    end do

    call halo_exchange_3d(VDnodes,dstruct) ! Qmax and Qmin could be synced here
    call halo_exchange_2d(Lengths,dstruct) ! Qmax and Qmin could be synced here

    do jnode=1,dstruct%nb_nodes
      VDnodes(:,jnode,:) = VDnodes(:,jnode,:) / Lengths(jnode)
    end do

    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      V(:,jedge,:) = 0.5_jprw * (VDnodes(:,ip1,:) + VDnodes(:,ip2,:) )
    end do

    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      V(:,iedge,YY) = 0.
    enddo

    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limited .and. (order .ge. 2) ) then
      Qmax(:,:) = -1e10
      Qmin(:,:) =  1e10
      call compute_Qmax_and_Qmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,Sx,Sy,Vx,Vy,ip1,ip2,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      do jlev=1,dstruct%nb_levels
        aun(jlev,jedge) = VDS(jlev,jedge)
        apos = max(0._jprw,aun(jlev,jedge))
        aneg = min(0._jprw,aun(jlev,jedge))
        fluxv(jlev,jedge) = Q(jlev,ip1)*apos + Q(jlev,ip2)*aneg
      end do
    enddo
    !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          do jlev=1,dstruct%nb_levels
            adv = adv + dstruct%sign(jedge,jnode)*fluxv(jlev,iedge)
          end do
        enddo
      endif
      ! Update the unknowns in vertices
      do jlev=1,dstruct%nb_levels
        Q(jlev,jnode) = Q(jlev,jnode) - adv/max(eps,volD(jlev,jnode)) * dt
        Qtmp(jlev,jnode) = Q(jlev,jnode)
        Q(jlev,jnode) = Q(jlev,jnode)*DR(jlev,jnode)
      end do
    enddo
    !$OMP END PARALLEL DO

    call halo_exchange_3d(Qtmp,dstruct) ! Qmax and Qmin could be synced here
    call halo_exchange_3d(Q,dstruct) ! Qmax and Qmin could be synced here





    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector_component, dstruct)

      call halo_exchange_3d(gradQ,dstruct)

      ! Compute antidiffusive normal velocity in faces

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dQdx,dQdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        do jlev=1,dstruct%nb_levels
          ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
          volume_of_two_cells = max(eps, volD(jlev,ip1) + volD(jlev,ip2) )
          dQdx = (gradQ(jlev,ip1,XX)+gradQ(jlev,ip2,XX)) / volume_of_two_cells
          dQdy = (gradQ(jlev,ip1,YY)+gradQ(jlev,ip2,YY)) / volume_of_two_cells
          Vx = V(jlev,jedge,XX)
          Vy = V(jlev,jedge,YY)
          ! variable sign option with asymptotic analysis, (mpdata gauge)
          aun(jlev,jedge) = abs(aun(jlev,jedge))*(Q(jlev,ip2)-Q(jlev,ip1))*0.5_jprw &
            &          -0.5_jprw*dt*aun(jlev,jedge)*(Vx*dQdx+Vy*dQdy)  ! = VDS*dQds
        end do
      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limited) then
        !call compute_Qmax_and_Qmin()
        call limit_antidiffusive_velocity()
      endif

      ! Compute fluxes from (limited) antidiffusive velocity
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
      do jnode=1,dstruct%nb_nodes
        adv = 0.0
        if(dstruct%nb_neighbours(jnode) > 1) then
          do jedge = 1,dstruct%nb_neighbours(jnode)
            iedge = dstruct%my_edges(jedge,jnode)
            do jlev=1,dstruct%nb_levels
              adv = adv + dstruct%sign(jedge,jnode)*aun(jlev,iedge)
            end do
          enddo
        endif
        ! Update the unknowns in vertices
        do jlev=1,dstruct%nb_levels
          Q(jlev,jnode) = Qtmp(jlev,jnode) - adv/max(eps,volD(jlev,jnode)) * dt
          Q(jlev,jnode) = Q(jlev,jnode) * DR(jlev,jnode)
        end do
      enddo
      !$OMP END PARALLEL DO

      call halo_exchange_3d(Q,dstruct)

    end do ! other passes

  contains
    
    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1, Q2    
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,Q1,Q2)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels
          Q1 = Q(jlev,jnode)
          Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q1 )
          Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q1 )
        end do
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          do jlev=1,dstruct%nb_levels
            Q2 = Q(jlev,dstruct%neighbours(jedge,jnode))
            if (.not. Q_is_vector_component) then
              Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q2 )
              Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q2 )
            else
              Qmax(jlev,jnode) = max( Qmax(jlev,jnode), pole_bc(iedge)*Q2 )
              Qmin(jlev,jnode) = min( Qmin(jlev,jnode), pole_bc(iedge)*Q2 )
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(Qmin,dstruct)
      call halo_exchange_3d(Qmax,dstruct)

    end subroutine compute_Qmax_and_Qmin

    subroutine limit_antidiffusive_velocity

      real(kind=jprw) :: asignp,asignn
      real(kind=jprw) :: limit = 1.  ! 1: second order, 0: first order

      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,rhin,rhout,jedge,iedge,apos,aneg,asignp,asignn)
      do jnode=1,dstruct%nb_nodes
        rhin  = 0.
        rhout = 0.
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          do jlev=1,dstruct%nb_levels
            apos = max(0._jprw,aun(jlev,iedge))
            aneg = min(0._jprw,aun(jlev,iedge))
            rhin(jlev)  = rhin(jlev)  - asignp*aneg - asignn*apos
            rhout(jlev) = rhout(jlev) + asignp*apos + asignn*aneg
          end do
        end do
        do jlev=1,dstruct%nb_levels
          cp(jlev,jnode) = ( Qmax(jlev,jnode)-Q(jlev,jnode) )*volD(jlev,jnode)/( rhin(jlev) * dt + eps )
          cn(jlev,jnode) = ( Q(jlev,jnode)-Qmin(jlev,jnode) )*volD(jlev,jnode)/( rhout(jlev)* dt + eps )
        end do
      end do
      !$OMP END PARALLEL DO

      call halo_exchange_3d(cp,dstruct)
      call halo_exchange_3d(cn,dstruct)

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)
        do jlev=1,dstruct%nb_levels
          if(aun(jlev,jedge) > 0._jprw) then
            aun(jlev,jedge)=aun(jlev,jedge)*min(limit,cp(jlev,ip2),cn(jlev,ip1))
          else
            aun(jlev,jedge)=aun(jlev,jedge)*min(limit,cn(jlev,ip2),cp(jlev,ip1))
          end if
        end do
      end do
      !$OMP END PARALLEL DO

    end subroutine limit_antidiffusive_velocity

  end subroutine mpdata_gauge_Q
  
  
  
  
  

  subroutine compute_gradient(Q,gradQ,Q_is_vector_component,dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(in)    :: Q(:,:)
    real(kind=jprw), intent(inout) :: gradQ(:,:,:)
    logical, intent(in) :: Q_is_vector_component
    real(kind=jprw), pointer :: S(:,:)
    real(kind=jprw) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2,jnode,jlev
    real(kind=jprw) :: avgQSx(dstruct%nb_levels,dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(dstruct%nb_levels,dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)

    ! derivatives 

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(XX,jedge)
      Sy  = S(YY,jedge)
      do jlev=1,dstruct%nb_levels
        avgQSx(jlev,jedge) = Sx*( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw
        avgQSy(jlev,jedge) = Sy*( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw
      end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(:,:,jnode) = 0.
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        do jlev=1,dstruct%nb_levels
          gradQ(XX,jlev,jnode) = gradQ(XX,jlev,jnode)+dstruct%sign(jedge,jnode)*avgQSx(jlev,iedge)
          gradQ(YY,jlev,jnode) = gradQ(YY,jlev,jnode)+dstruct%sign(jedge,jnode)*avgQSy(jlev,iedge)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! special treatment for the north & south pole cell faces
    ! Sx == 0 at pole, and Sy has same sign at both sides of pole
    if (.not. Q_is_vector_component) then
      do jedge = 1,dstruct%nb_pole_edges
        iedge = dstruct%pole_edges(jedge)
        ip1   = dstruct%edges(iedge,1)
        ip2   = dstruct%edges(iedge,2)
        Sy    = S(YY,iedge)

        do jlev=1,dstruct%nb_levels
          avgQ  = ( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw

          ! correct for wrong Y-derivatives in previous loop,
          gradQ(YY,jlev,ip2) = gradQ(YY,jlev,ip2) + 2._jprw*Sy*avgQ
        end do
      end do
    end if
  end subroutine compute_gradient

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
  use mpdata3d_module, &
    & only: mpdata_gauge_D, mpdata_gauge_Q, &
    &       compute_gradient

  implicit none
  private
  public :: setup_isentropic
  public :: propagate_state
  public :: set_state_zonal_flow
  public :: set_topography
  public :: set_time_step

  real(kind=jprw), parameter :: eps    = 1.e-6
  real(kind=jprw), parameter :: radius = 63.6620e+03 !6371.22e+03
  real(kind=jprw), parameter :: f0     = 1.4584e-04 !coriolis parameter (=2xearth's omega)
  real(kind=jprw), parameter :: grav   = 9.80616
  real(kind=jprw), parameter :: pi     = acos(-1._jprw)

  real(kind=jprw), parameter :: Rgas   = 287.04 ! gas constant
  real(kind=jprw), parameter :: stf    = 1.02e-05 !stf=bv**2/g
  real(kind=jprw), parameter :: cp     = 3.5*Rgas
  real(kind=jprw), parameter :: cap    = Rgas/cp
  real(kind=jprw), parameter :: pscal  = 1.e5
  real(kind=jprw), parameter :: pscali = 1.e-5
  real(kind=jprw), parameter :: dz     = 35.
  real(kind=jprw), parameter :: pr00   = 1.e5
  real(kind=jprw), parameter :: th00   = 293.15

  real(kind=jprw), parameter :: ibs=0.     ! 0: isentropic/ isosteric --- 1: isopicnic
  real(kind=jprw), parameter :: icp=1.-ibs

  real(kind=jprw), allocatable :: alfs(:), alf0(:), z0(:)


  real(kind=jprw) :: dt_forward = 20.
  integer :: iter = 0



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
    real(kind=jprw), pointer :: coords(:,:), vol(:,:), hx(:), hy(:), dhxdy_over_G(:,:), pole_bc(:), G(:,:)
    integer :: jnode, jedge, iedge, jlev

    call create_field_in_nodes_3d("jacobian",1,dstruct)
    call create_field_in_nodes_2d("hx",1,dstruct)
    call create_field_in_nodes_2d("hy",1,dstruct)
    call create_field_in_nodes_3d("dhxdy_over_G",1,dstruct)
    call create_field_in_edges_2d("pole_bc",1,dstruct)


    coords => vector_field_2d("coordinates",dstruct)
    vol    => scalar_field_3d("dual_volumes",dstruct)
    hx    => scalar_field_2d("hx",dstruct)
    hy    => scalar_field_2d("hy",dstruct)
    G     => scalar_field_3d("jacobian",dstruct)
    dhxdy_over_G => scalar_field_3d("dhxdy_over_G",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)
    !dir$ ivdep
    do jnode=1,dstruct%nb_nodes
      y = coords(YY,jnode)
      cos_y = cos(y)
      sin_y = sin(y)
      hx(jnode) = radius*cos_y
      hy(jnode) = radius
      do jlev=1,dstruct%nb_levels
        G(jlev,jnode) = hx(jnode)*hy(jnode)
        vol(jlev,jnode) = vol(jlev,jnode)*G(jlev,jnode)
        dhxdy_over_G(jlev,jnode) = - sin_y/(radius*max(eps,cos_y))
      end do
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

    call log_debug("Backup solution")
    call backup_solution(dstruct)

    if (step == 0) then ! Pre-compute forcing
      call log_debug("Compute forcing")
      call compute_forcing(dstruct)

      call log_debug("Compute advective velocities")
      call compute_advective_velocities(dt,dstruct,"extrapolate")
    end if
    
    call log_debug("Add fording to solution")
    call add_forcing_to_solution(dt,dstruct)
    
    call log_debug("Advect solution")
    call advect_solution(dt,dstruct)

    call log_debug("Implicit solve")
    call implicit_solve(dt,dstruct)

    call log_debug("Compute advective velocities")
    call compute_advective_velocities(dt,dstruct,"extrapolate")

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

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
    do jnode=1,dstruct%nb_nodes
      x=coords(XX,jnode)
      y=coords(YY,jnode)
      cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
      U0(XX,jnode) =  omega*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*max(0._jprw,cos(y))
      U0(YY,jnode) = -omega*sin(x)*sin(beta)*radius

      dlh(jnode)   = -radius**2*(f0+omega)*0.5_jprw*omega/grav &
                   & *(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2

      ! Initial Z at levels
      do jlev=1,dstruct%nb_levels
        Z0(jlev) = (jlev-1)*dz
        alf0(jlev) = alfb(Z0(jlev),0._jprw,stf)
      end do
      ! Z at half levels (staggered)
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
        D_amb(jlev,jnode)    = (press(jlev,jnode)-press(jlev+1,jnode))!*pscali
        Q_amb(XX,jlev,jnode) = (U0(XX,jnode)-0.01*cos(y)*H(jlev,jnode))*D_amb(jlev,jnode)
        Q_amb(YY,jlev,jnode) = U0(YY,jnode)*D_amb(jlev,jnode)
        D(jlev,jnode)   = D_amb(jlev,jnode)
        Q(:,jlev,jnode) = Q_amb(:,jlev,jnode)
        U(:,jlev,jnode) = Q(:,jlev,jnode)/max(eps,D(jlev,jnode))
      enddo
    end do
    !$OMP END PARALLEL DO


    call compute_forcing(dstruct)

  end subroutine set_state_zonal_flow


  subroutine set_topography(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:), pointer :: H0
    real(kind=jprw), dimension(:,:), pointer :: coords
    real(kind=jprw) :: zlatc, zlonc, zr, zz, rad, x, y, amp
    integer :: jnode

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
  end subroutine set_topography


  subroutine backup_solution(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: D(:,:), D0(:,:)
    real(kind=jprw), pointer :: Q(:,:,:), Q0(:,:,:)
    integer :: jnode
    
    D  => scalar_field_3d("depth",dstruct)
    D0 => scalar_field_3d("depth_backup",dstruct)
    Q  => vector_field_3d("velocity",dstruct)
    Q0 => vector_field_3d("velocity_backup",dstruct)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      D0(:,jnode)   = D(:,jnode)
      Q0(:,jnode,:) = Q(:,jnode,:)
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
    real(kind=jprw) :: Ux, Uy, U0x, U0y, Vx, Vy, Rx, Ry, dVxdx, dVxdy, dVydx, dVydy
    integer :: jnode, jedge, iedge, jlev, ip1, ip2
    real(kind=jprw), dimension(:),   pointer :: hx, hy
    real(kind=jprw), dimension(:,:),   pointer :: D, D0, vol, coords
    real(kind=jprw), dimension(:,:,:), pointer :: U, U0, R, Vedges
    real(kind=jprw) :: Vnodes(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw) :: grad_Vx(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw) :: grad_Vy(dstruct%nb_levels,dstruct%nb_nodes,2)

    coords => vector_field_2d("coordinates",dstruct)
    Vedges => vector_field_3d("advective_velocity",dstruct)
    D      => scalar_field_3d("depth",dstruct)
    D0     => scalar_field_3d("depth_backup",dstruct)
    U      => vector_field_3d("velocity",dstruct)
    U0     => vector_field_3d("velocity_backup",dstruct)
    hx     => scalar_field_2d("hx",dstruct)
    hy     => scalar_field_2d("hy",dstruct)
    R      => vector_field_3d("velocity_forcing",dstruct)
    vol    => scalar_field_3d("dual_volumes",dstruct)

    if( option .eq. "advect") then
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
      do jnode=1,dstruct%nb_nodes
        Vnodes(:,jnode,XX)=(U(:,jnode,XX)+0.5*dt*R(:,jnode,XX)) 
        Vnodes(:,jnode,YY)=(U(:,jnode,YY)+0.5*dt*R(:,jnode,YY))
      end do
      !$OMP END PARALLEL DO
      call compute_gradient( Vnodes(:,:,XX), grad_Vx, .True., dstruct )
      call compute_gradient( Vnodes(:,:,YY), grad_Vy, .True., dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Vx,Vy,Rx,Ry,dVxdx,dVxdy,dVydx,dVydy)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels
          Vx = Vnodes(jlev,jnode,XX)
          Vy = Vnodes(jlev,jnode,YY)
          Rx = R(jlev,jnode,XX)
          Ry = R(jlev,jnode,YY)
          dVxdx = grad_Vx(jlev,jnode,XX)*hy(jnode)/vol(jlev,jnode)
          dVxdy = grad_Vx(jlev,jnode,YY)*hx(jnode)/vol(jlev,jnode)
          dVydx = grad_Vy(jlev,jnode,YY)*hy(jnode)/vol(jlev,jnode)    
          dVydy = grad_Vy(jlev,jnode,YY)*hx(jnode)/vol(jlev,jnode)

          Vnodes(jlev,jnode,XX) = ( Vx - 0.5*dt*(Vx*dVxdx+Vy*dVxdy)) * hy(jnode)
          Vnodes(jlev,jnode,YY) = ( Vy - 0.5*dt*(Vx*dVydx+Vy*dVydy)) * hx(jnode)
        end do
      end do
      !$OMP END PARALLEL DO
      call halo_exchange_3d( Vnodes, dstruct )
    else if( option .eq. "extrapolate") then
      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy,U0x,U0y)
      do jnode=1,dstruct%nb_nodes
        do jlev=1,dstruct%nb_levels
          Ux    = U(jlev,jnode,XX)
          Uy    = U(jlev,jnode,YY)
          U0x   = U0(jlev,jnode,XX)
          U0y   = U0(jlev,jnode,YY)
          Vnodes(jlev,jnode,XX) = ( 1.5_jprw*Ux - 0.5_jprw*U0x ) * hy(jnode)
          Vnodes(jlev,jnode,YY) = ( 1.5_jprw*Uy - 0.5_jprw*U0y ) * hx(jnode)
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Vedges(:,jedge,XX) = (Vnodes(:,ip1,XX)+Vnodes(:,ip2,XX))*0.5_jprw
      Vedges(:,jedge,YY) = (Vnodes(:,ip1,YY)+Vnodes(:,ip2,YY))*0.5_jprw
    enddo
    !$OMP END PARALLEL DO
    
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      Vedges(:,iedge,YY) = 0.
    enddo

  end subroutine compute_advective_velocities

  subroutine compute_forcing(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, jlev
    real(kind=jprw) :: Qx, Qy
    real(kind=jprw), dimension(:),   pointer :: H0, hx, hy, cor, p0
    real(kind=jprw), dimension(:,:), pointer :: D, vol, dhxdy_over_G, coords, M, press
    real(kind=jprw), dimension(:,:,:), pointer :: Q, R
    real(kind=jprw) :: grad_M(2,dstruct%nb_levels,dstruct%nb_nodes)

    coords => vector_field_2d("coordinates",dstruct)
    vol => scalar_field_3d("dual_volumes",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_3d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    H0 => scalar_field_2d("topography",dstruct)
    D  => scalar_field_3d("depth",dstruct)
    Q  => vector_field_3d("momentum",dstruct)
    R  => vector_field_3d("forcing",dstruct)
    M  => scalar_field_3d("montgomery_potential",dstruct)
    press  => scalar_field_3d("pressure",dstruct)
    p0 => scalar_field_2d("p0",dstruct)

    do jnode=1,dstruct%nb_nodes
      ! Convert pressure to exner in case ibs==0 (isentropic)
      press(dstruct%nb_levels,jnode) = ibs*p0(jnode)+icp*cp*(p0(jnode)/pr00)**cap

      ! Define Montegomery potential (alias Bernoulli head)
      M(1,jnode) = alf0(1)*press(1,jnode)+grav*H0(jnode) + (alfs(1)-alf0(1))*press(1,jnode)

      do jlev=2,dstruct%nb_levels-1
         M(jlev,jnode) = M(jlev-1,jnode) + (alfs(jlev)-alf0(jlev-1))*press(jlev,jnode)
      end do
    end do

    call compute_gradient( M, grad_M, .False., dstruct )

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
    do jnode=1,dstruct%nb_nodes
      do jlev=1,dstruct%nb_levels-1
        Qx    = Q(XX,jlev,jnode)
        Qy    = Q(YY,jlev,jnode)
        R(XX,jlev,jnode) = -grav*D(jlev,jnode)*grad_M(XX,jlev,jnode)*hy(jnode)/vol(jlev,jnode) &
          &              + cor(jnode)*Qy - dhxdy_over_G(jlev,jnode)*Qx*Qy/max(eps,D(jlev,jnode))
        R(YY,jlev,jnode) = -grav*D(jlev,jnode)*grad_M(YY,jlev,jnode)*hx(jnode)/vol(jlev,jnode) &
          &              - cor(jnode)*Qx + dhxdy_over_G(jlev,jnode)*Qx*Qx/max(eps,D(jlev,jnode))
      end do
    end do
    !$OMP END PARALLEL DO

    call halo_exchange_3d(R,dstruct)

  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: U(:,:,:) , R(:,:,:), D(:,:)
    integer :: jnode

    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    R => vector_field_3d("velocity_forcing",dstruct)
    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      U(:,jnode,XX) = ( U(:,jnode,XX) + 0.5_jprw*dt*R(:,jnode,XX) )
      U(:,jnode,YY) = ( U(:,jnode,YY) + 0.5_jprw*dt*R(:,jnode,YY) )  
    end do
    !$OMP END PARALLEL DO
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, jlev, m
    real(kind=jprw) :: Ux, Uy, Rx, Ry
    real(kind=jprw) :: Ux_adv, Uy_adv, Rx_exp, Ry_exp

    real(kind=jprw), dimension(:),     pointer :: H0, hx, hy, cor
    real(kind=jprw), dimension(:,:),   pointer :: H, D, vol, dhxdy_over_G
    real(kind=jprw), dimension(:,:),   pointer :: coords
    real(kind=jprw), dimension(:,:,:), pointer :: U, R
    real(kind=jprw) :: grad_H(dstruct%nb_levels,dstruct%nb_nodes, 2)

    coords => vector_field_2d("coordinates",dstruct)
    vol => scalar_field_3d("dual_volumes",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_3d("height",dstruct)
    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    R => vector_field_3d("velocity_forcing",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_3d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    do jlev=1,dstruct%nb_levels
      H(jlev,:) = H0(:) + D(jlev,:)
    end do

    call compute_gradient( H, grad_H, .False., dstruct )

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy,Rx_exp,Ry_exp,Ux_adv,Uy_adv,m,Rx,Ry)
    do jnode=1,dstruct%nb_nodes
      do jlev=1,dstruct%nb_levels
        Ux    = U(jlev,jnode,XX)
        Uy    = U(jlev,jnode,YY)

        Rx_exp = -grav*grad_H(jlev,jnode,XX)*hy(jnode)/vol(jlev,jnode)
        Ry_exp = -grav*grad_H(jlev,jnode,YY)*hx(jnode)/vol(jlev,jnode)

        Ux_adv = Ux
        Uy_adv = Uy

        do m=1,3 ! Three iterations at most is enough to converge
          Rx = Rx_exp + cor(jnode)*Uy - dhxdy_over_G(jlev,jnode)*Ux*Uy
          Ry = Ry_exp - cor(jnode)*Ux + dhxdy_over_G(jlev,jnode)*Ux*Ux
          Ux = Ux_adv + 0.5_jprw*dt*Rx
          Uy = Uy_adv + 0.5_jprw*dt*Ry
        end do
        U(jlev,jnode,XX) = Ux
        U(jlev,jnode,YY) = Uy
        R(jlev,jnode,XX) = Rx_exp + cor(jnode)*Uy - dhxdy_over_G(jlev,jnode)*Ux*Uy
        R(jlev,jnode,YY) = Ry_exp - cor(jnode)*Ux + dhxdy_over_G(jlev,jnode)*Ux*Ux
      end do
    end do
    !$OMP END PARALLEL DO
    call halo_exchange_3d(U,dstruct)
    call halo_exchange_3d(R,dstruct)

  end subroutine implicit_solve



  subroutine advect_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:,:),   pointer :: D, D0, DR
    real(kind=jprw), dimension(:,:,:), pointer :: U, V
    real(kind=jprw) :: VDS(dstruct%nb_levels,dstruct%nb_edges)
    integer :: jnode
    
    D => scalar_field_3d("depth",dstruct)
    D0 => scalar_field_3d("depth_backup",dstruct)
    U => vector_field_3d("velocity",dstruct)
    V => vector_field_3d("advective_velocity",dstruct)
    DR => scalar_field_3d("depth_ratio",dstruct)
   
    !    mpdata_gauge_D( time, variable, velocity, VDS,  order, limit,  is_vector, dstruct )
    call mpdata_gauge_D( dt,   D,        V,        VDS,  1,     .True., .False.,   dstruct )
    
    ! compute ratio
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      DR(:,jnode) = D0(:,jnode) / max( D(:,jnode), eps )
    end do
    !$OMP END PARALLEL DO

    !    mpdata_gauge_Q( time, variable, VDS, DR, D0,  order, limit,  is_vector, dstruct )
    call mpdata_gauge_Q( dt,   U(:,:,XX),  VDS, DR, D0,  1,     .True., .True. ,   dstruct )
    call mpdata_gauge_Q( dt,   U(:,:,YY),  VDS, DR, D0,  1,     .True., .True. ,   dstruct )
  end subroutine advect_solution

end module isentropic_module
