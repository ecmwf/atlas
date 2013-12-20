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
  public :: mpdata_gauge_D
  public :: mpdata_gauge_U
  public :: compute_gradient
  public :: compute_gradient_tensor

  real(kind=jprw), parameter :: eps = 1.e-6

contains

  subroutine mpdata_gauge_D(dt,D,V,VDS,order,limited,dstruct)
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: D(:)
    real(kind=jprw), intent(in) :: V(:,:)
    logical, intent(in) :: limited
    integer, intent(in) :: order
    real(kind=jprw), intent(out) :: VDS(dstruct%nb_edges)

    integer :: jnode, jedge, iedge, jpass, ip1,ip2
    real(kind=jprw) :: sx, sy, volume_of_two_cells, dDdx, dDdy, Vx, Vy, apos, aneg
    real(kind=jprw) :: Dmin(dstruct%nb_nodes)
    real(kind=jprw) :: Dmax(dstruct%nb_nodes)
    real(kind=jprw) :: rhin
    real(kind=jprw) :: rhout
    real(kind=jprw) :: cp(dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(dstruct%nb_edges)
    real(kind=jprw) :: gradD(2,dstruct%nb_nodes)
    real(kind=jprw), pointer :: vol(:), S(:,:), pole_bc(:)

    vol     => scalar_field_2d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)


    VDS(:) = 0.

    ! 1. First pass
    ! -------------
    jpass = 1

    ! non-oscillatory option
    if( limited .and. (order .ge. 2) ) then
      Dmax(:) = -1e10
      Dmin(:) =  1e10
      call compute_Dmax_and_Dmin()
    end if

    ! Compute the normal velocity in faces, and advection in vertices

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,Sx,Sy,Vx,Vy,ip1,ip2,apos,aneg)
    do jedge = 1,dstruct%nb_edges
      Sx = S(XX,jedge)
      Sy = S(YY,jedge)
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

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
     D(jnode) = D(jnode) - adv/vol(jnode) * dt
    enddo
    !$OMP END PARALLEL DO

    call halo_exchange(D,dstruct) ! Dmax and Dmin could be synced here

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient(D, gradD, dstruct)

      call halo_exchange(gradD,dstruct)

      ! Compute antidiffusive normal velocity in faces

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
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(jedge) = abs(aun(jedge))*(D(ip2)-D(ip1))*0.5_jprw &
          &          -0.5_jprw*dt*aun(jedge)*(Vx*dDdx+Vy*dDdy)
      end do
      !$OMP END PARALLEL DO

      ! non-oscillatory option
      if (limited) then
        call limit_antidiffusive_velocity()
      endif

      do jedge = 1,dstruct%nb_edges
        VDS(jedge) = VDS(jedge) + aun(jedge)
      end do

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
        D(jnode) = max( D(jnode) - adv/vol(jnode) * dt, 0. )
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
      real(kind=jprw) :: limit = 1.  ! 1: second order, 0: first order

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

  end subroutine mpdata_gauge_D
  
  
  
  
  
  
  
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

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
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


  subroutine compute_gradient_tensor(Q,gradQ,dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(in)  :: Q(:,:)
    real(kind=jprw), intent(out) :: gradQ(:,:)
    real(kind=jprw), pointer :: S(:,:)
    real(kind=jprw) :: Sx,Sy
    integer :: jedge,iedge,ip1,ip2,jnode
    real(kind=jprw) :: avgQSx(2,dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(2,dstruct%nb_edges)

    S   => vector_field_2d("dual_normals",dstruct)

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
    & only: mpdata_gauge_D, mpdata_gauge_U, &
    &       compute_gradient, compute_gradient_tensor

  implicit none
  private
  public :: setup_shallow_water 
  public :: propagate_state
  public :: set_state_rossby_haurwitz
  public :: set_state_zonal_flow
  public :: set_topography
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

  subroutine setup_shallow_water(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw) :: y, cos_y, sin_y
    real(kind=jprw), pointer :: coords(:,:), vol(:), hx(:), hy(:), dhxdy_over_G(:), pole_bc(:), G(:)
    integer :: jnode, jedge, iedge
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
    call create_field_in_nodes_2d("velocity_forcing",2,dstruct)
    call create_field_in_nodes_2d("depth_backup",1,dstruct)
    call create_field_in_nodes_2d("velocity_backup",2,dstruct)

    call create_field_in_edges_2d("advective_velocity",2,dstruct)
    call create_field_in_nodes_2d("depth_ratio",1,dstruct) ! old/new
    
    call create_field_in_nodes_2d("topography",1,dstruct)
    call create_field_in_nodes_2d("height",1,dstruct)

  end subroutine setup_shallow_water

  subroutine propagate_state(dt,order,dstruct)
    real(kind=jprw), intent(in) :: dt
    integer, intent(in) :: order
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
      call step_forward(iter,dt_fwd,order,dstruct)

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

  end subroutine propagate_state




  subroutine step_forward(step,dt,order,dstruct)
    integer, intent(inout) :: step
    real(kind=jprw), intent(in) :: dt
    integer, intent(in) :: order
    type(DataStructure_type), intent(inout) :: dstruct

    call backup_solution(dstruct)

    if (step == 0) then ! Pre-compute forcing
      call compute_forcing(dstruct)

      call compute_advective_velocities(dt,dstruct,"extrapolate")
    end if
    
    call add_forcing_to_solution(dt,dstruct)
    
    call advect_solution(dt,order,dstruct)

    call implicit_solve(dt,dstruct)

    call compute_advective_velocities(dt,dstruct,"extrapolate")

    dstruct%time = dstruct%time+dt
    step = step+1

  end subroutine step_forward




  subroutine set_state_rossby_haurwitz(dstruct)
    type(DataStructure_type), intent(inout)      :: dstruct
    real(kind=jprw), dimension(:), pointer   :: D, cor, H0, H
    real(kind=jprw), dimension(:,:), pointer :: U, coords
    integer :: jnode, ir
    real(kind=jprw) :: aaa0,zk,om,ph0,x,y, sin_y, cos_y

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3
    aaa0 = 0.

    coords => vector_field_2d("coordinates",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
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
      D(jnode) = (ph0+radius**2*fa(y)+radius**2*fb(y)*cos(ir*x) &
        &        +radius**2*fc(y)*cos(2._jprw*ir*x)) / grav
      D(jnode) = max(aaa0,D(jnode) - H0(jnode))
      if(y == 0.5_jprw*pi) U(XX,jnode)=0.
      if(y ==-0.5_jprw*pi) U(XX,jnode)=0.
      H(jnode) = H0(jnode) + D(jnode)
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
    type(DataStructure_type), intent(inout)      :: dstruct
    real(kind=jprw), dimension(:), pointer   :: D, cor, H0, H
    real(kind=jprw), dimension(:,:), pointer :: U, coords
    integer :: jnode
    real(kind=jprw) :: x,y
    real(kind=jprw), parameter :: USCAL = 20.
    real(kind=jprw), parameter :: H00 = grav * 8e3
    real(kind=jprw), parameter :: pvel = USCAL/radius
    real(kind=jprw), parameter :: beta = 0.! pi/4._jprw


    coords => vector_field_2d("coordinates",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_2d("height",dstruct)

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
    do jnode=1,dstruct%nb_nodes
      x=coords(XX,jnode)
      y=coords(YY,jnode)
      cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
      D(jnode)     = (H00-radius**2*(f0+pvel)*0.5*pvel*(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2)
      D(jnode)     = max(0._jprw, D(jnode)/grav - H0(jnode))
      U(XX,jnode)  =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y)
      U(YY,jnode)  = -pvel*sin(x)*sin(beta)*radius
      H(jnode) = H0(jnode) + D(jnode)
    end do
    !$OMP END PARALLEL DO

  end subroutine set_state_zonal_flow


  subroutine set_topography(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:), pointer :: H0
    real(kind=jprw), dimension(:,:), pointer :: coords
    real(kind=jprw) :: amp = 0. ! 2000. ! amplitude of hill
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
      if (dist.le.rad) then
        H0(jnode) = amp * (1.-gamm*dist/rad)
      else
        H0(jnode) = 0.
      end if
    end do
  end subroutine set_topography


  subroutine backup_solution(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:),   pointer :: D, D0
    real(kind=jprw), dimension(:,:), pointer :: Q, Q0
    integer :: jnode
    
    D  => scalar_field_2d("depth",dstruct)
    D0 => scalar_field_2d("depth_backup",dstruct)
    Q  => vector_field_2d("velocity",dstruct)
    Q0 => vector_field_2d("velocity_backup",dstruct)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      D0(jnode)   = D(jnode)
      Q0(:,jnode) = Q(:,jnode)
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
    integer :: jnode, jedge, iedge, ip1, ip2
    real(kind=jprw), dimension(:),   pointer :: D, D0, hx, hy, vol
    real(kind=jprw), dimension(:,:), pointer :: U, U0, R, Vedges, coords
    real(kind=jprw) :: Vnodes(2,dstruct%nb_nodes), grad_Vnodes(4,dstruct%nb_nodes)

    coords => vector_field_2d("coordinates",dstruct)
    Vedges => vector_field_2d("advective_velocity",dstruct)
    D      => scalar_field_2d("depth",dstruct)
    D0     => scalar_field_2d("depth_backup",dstruct)
    U      => vector_field_2d("velocity",dstruct)
    U0     => vector_field_2d("velocity_backup",dstruct)
    hx     => scalar_field_2d("hx",dstruct)
    hy     => scalar_field_2d("hy",dstruct)
    R      => vector_field_2d("velocity_forcing",dstruct)
    vol    => scalar_field_2d("dual_volumes",dstruct)

    if( option .eq. "advect") then
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
      do jnode=1,dstruct%nb_nodes
        Vnodes(:,jnode)=(U(:,jnode)+0.5*dt*R(:,jnode)) 
      end do
      !$OMP END PARALLEL DO
      call compute_gradient_tensor( Vnodes, grad_Vnodes, dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Vx,Vy,Rx,Ry,dVxdx,dVxdy,dVydx,dVydy)
      do jnode=1,dstruct%nb_nodes
        Vx = Vnodes(XX,jnode)
        Vy = Vnodes(YY,jnode)
        Rx = R(XX,jnode)
        Ry = R(YY,jnode)
        dVxdx = grad_Vnodes(XXDXX,jnode)*hy(jnode)/vol(jnode)
        dVxdy = grad_Vnodes(XXDYY,jnode)*hx(jnode)/vol(jnode)
        dVydx = grad_Vnodes(YYDXX,jnode)*hy(jnode)/vol(jnode)    
        dVydy = grad_Vnodes(YYDYY,jnode)*hx(jnode)/vol(jnode)

        Vnodes(XX,jnode) = ( Vx - 0.5*dt*(Vx*dVxdx+Vy*dVxdy)) * hy(jnode)
        Vnodes(YY,jnode) = ( Vy - 0.5*dt*(Vx*dVydx+Vy*dVydy)) * hx(jnode)
      enddo
      !$OMP END PARALLEL DO
      call halo_exchange( Vnodes, dstruct )
    else if( option .eq. "extrapolate") then
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
    end if

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      VedgeS(:,jedge) = (Vnodes(:,ip1)+Vnodes(:,ip2))*0.5_jprw
    enddo
    !$OMP END PARALLEL DO
    
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      VedgeS(YY,iedge) = 0.
    enddo

  end subroutine compute_advective_velocities

  subroutine compute_forcing(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode
    real(kind=jprw) :: Ux, Uy
    real(kind=jprw), dimension(:),   pointer :: H, H0, D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:), pointer :: U, R
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
    R => vector_field_2d("velocity_forcing",dstruct)

    H(:) = H0(:) + D(:)
    call compute_gradient( H, grad_H, dstruct )

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

    call halo_exchange(R,dstruct)
  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: U(:,:) , R(:,:), D(:)
    integer :: jnode
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    R => vector_field_2d("velocity_forcing",dstruct)
    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      U(:,jnode) = ( U(:,jnode) + 0.5_jprw*dt*R(:,jnode) )
    end do
    !$OMP END PARALLEL DO
  end subroutine add_forcing_to_solution



  subroutine implicit_solve(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, m
    real(kind=jprw) :: Ux, Uy, Rx, Ry
    real(kind=jprw) :: Ux_adv, Uy_adv, Rx_exp, Ry_exp

    real(kind=jprw), dimension(:),   pointer :: H, H0, D, vol, hx, hy, dhxdy_over_G, cor
    real(kind=jprw), dimension(:,:), pointer :: U, R
    real(kind=jprw) :: grad_H(2,dstruct%nb_nodes)

    vol => scalar_field_2d("dual_volumes",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_2d("height",dstruct)
    D => scalar_field_2d("depth",dstruct)
    U => vector_field_2d("velocity",dstruct)
    R => vector_field_2d("velocity_forcing",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_2d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    ! D is already up to date at time level (n+1), just by MPDATA advection
    H(:) = H0(:) + D(:)
    call compute_gradient( H, grad_H, dstruct )

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
    end do
    !$OMP END PARALLEL DO
    call halo_exchange(U,dstruct)
    call halo_exchange(R,dstruct)

  end subroutine implicit_solve



  subroutine advect_solution(dt,order,dstruct)
    real(kind=jprw), intent(in) :: dt
    integer, intent(in) :: order
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:),   pointer :: D, D0, DR
    real(kind=jprw), dimension(:,:), pointer :: U, V
    real(kind=jprw) :: VDS(dstruct%nb_edges)
    integer :: jnode
    
    D => scalar_field_2d("depth",dstruct)
    D0 => scalar_field_2d("depth_backup",dstruct)
    U => vector_field_2d("velocity",dstruct)
    V => vector_field_2d("advective_velocity",dstruct)
    DR => scalar_field_2d("depth_ratio",dstruct)
   
    !    mpdata_gauge_D( time, variable, velocity, VDS,  order, limit,   dstruct )
    call mpdata_gauge_D( dt,   D,        V,        VDS,  order, .True.,  dstruct )

    ! compute ratio
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      DR(jnode) = D0(jnode) / max( D(jnode), eps )
    end do
    !$OMP END PARALLEL DO
    !    mpdata_gauge_U( time, variable, VDS, DR, D0,  order, limit,  dstruct )
    call mpdata_gauge_U( dt,   U,        VDS, DR, D0,  order, .True., dstruct )
  end subroutine advect_solution

end module shallow_water_module
