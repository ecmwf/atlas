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
    real(kind=jprw), intent(out) :: VDS(dstruct%nb_edges)
    
    integer :: jnode, jedge, iedge, jpass, ip1,ip2, jlev
    real(kind=jprw) :: sx, sy, volume_of_two_cells, dQdx, dQdy, Vx, Vy, apos, aneg
    real(kind=jprw) :: Qmin(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: Qmax(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: rhin
    real(kind=jprw) :: rhout
    real(kind=jprw) :: cp(dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(dstruct%nb_edges)
    real(kind=jprw) :: fluxv(dstruct%nb_edges)
    real(kind=jprw) :: gradQ(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw), pointer :: vol(:,:), S(:,:), pole_bc(:)
    

    vol     => scalar_field_3d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)

    jlev = 1

    VDS(:) = 0.

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
      Sx = S(jedge,XX)
      Sy = S(jedge,YY)
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)

      Vx = V(jlev,jedge,XX)
      Vy = V(jlev,jedge,YY)

      aun(jedge) = Vx*Sx + Vy*Sy 

      apos = max(0._jprw,aun(jedge))
      aneg = min(0._jprw,aun(jedge))
      fluxv(jedge) = Q(jlev,ip1)*apos + Q(jlev,ip2)*aneg
      VDS(jedge) = fluxv(jedge)
    enddo
    !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          adv = adv + dstruct%sign(jedge,jnode)*fluxv(iedge)
        enddo
      endif
     ! Update the unknowns in vertices
     Q(jlev,jnode) = Q(jlev,jnode) - adv/vol(jlev,jnode) * dt
    enddo
    !$OMP END PARALLEL DO

    call synchronise(Q,dstruct) ! Qmax and Qmin could be synced here

    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector_component, dstruct)

      call synchronise(gradQ,dstruct)

      ! Compute antidiffusive normal velocity in faces

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dQdx,dQdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = vol(jlev,ip1) + vol(jlev,ip2)
        dQdx = (gradQ(jlev,ip1,XX)+gradQ(jlev,ip2,XX)) / volume_of_two_cells
        dQdy = (gradQ(jlev,ip1,YY)+gradQ(jlev,ip2,YY)) / volume_of_two_cells
        Vx = V(jlev,jedge,XX)
        Vy = V(jlev,jedge,YY)
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(jedge) = abs(aun(jedge))*(Q(jlev,ip2)-Q(jlev,ip1))*0.5_jprw &
          &          -0.5_jprw*dt*aun(jedge)*(Vx*dQdx+Vy*dQdy)
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
            fluxv(iedge) = aun(iedge)
            adv = adv + dstruct%sign(jedge,jnode)*fluxv(iedge)
          enddo
        endif
        ! Update the unknowns in vertices
        Q(jlev,jnode) = max( Q(jlev,jnode) - adv/vol(jlev,jnode) * dt, 0. )
      enddo
      !$OMP END PARALLEL DO
      call synchronise(Q,dstruct)

    end do ! other passes

  contains
    
    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1, Q2    
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,Q1,Q2)
      do jnode=1,dstruct%nb_nodes
        Q1 = Q(jlev,jnode)
        Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q1 )
        Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q1 )
        do jedge = 1,dstruct%nb_neighbours(jnode)
          Q2 = Q(jlev,dstruct%neighbours(jedge,jnode))
          iedge = dstruct%my_edges(jedge,jnode)
          if (.not. Q_is_vector_component) then
            Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q2 )
            Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q2 )
          else
            Qmax(jlev,jnode) = max( Qmax(jlev,jnode), pole_bc(iedge)*Q2 )
            Qmin(jlev,jnode) = min( Qmin(jlev,jnode), pole_bc(iedge)*Q2 )
          end if
        end do
      end do
      !$OMP END PARALLEL DO

      call synchronise(Qmin,dstruct)
      call synchronise(Qmax,dstruct)

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
          apos = max(0._jprw,aun(iedge))
          aneg = min(0._jprw,aun(iedge))
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          rhin  = rhin  - asignp*aneg - asignn*apos
          rhout = rhout + asignp*apos + asignn*aneg
        end do
        cp(jnode) = ( Qmax(jlev,jnode)-Q(jlev,jnode) )*vol(jlev,jnode)/( rhin * dt + eps )
        cn(jnode) = ( Q(jlev,jnode)-Qmin(jlev,jnode) )*vol(jlev,jnode)/( rhout* dt + eps )
      end do
      !$OMP END PARALLEL DO

      call synchronise(cp,dstruct)
      call synchronise(cn,dstruct)

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
  
  
  
  
  
  
  
  subroutine mpdata_gauge_Q(dt,Q,VDS,DR,D,order,limited,Q_is_vector_component,dstruct)
    real(kind=jprw), intent(in)  :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(inout) :: Q(:,:), D(:,:)
    real(kind=jprw), intent(in) :: VDS(:), DR(:,:)
    logical, intent(in) :: Q_is_vector_component, limited
    integer, intent(in) :: order
    integer :: jnode, jedge, iedge, jpass, jlev, ip1,ip2
    real(kind=jprw) :: Sx, Sy, Ssqr, volume_of_two_cells, dQdx, dQdy, Vx, Vy
    real(kind=jprw) :: apos, aneg, x1, x2, y1, y2, length
    real(kind=jprw) :: Qmin(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: Qmax(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: rhin
    real(kind=jprw) :: rhout
    real(kind=jprw) :: cp(dstruct%nb_nodes)
    real(kind=jprw) :: cn(dstruct%nb_nodes)
    real(kind=jprw) :: adv
    real(kind=jprw) :: aun(dstruct%nb_edges)
    real(kind=jprw) :: fluxv(dstruct%nb_edges)
    real(kind=jprw) :: gradQ(dstruct%nb_levels,dstruct%nb_nodes,2)
    real(kind=jprw) :: V(dstruct%nb_edges,2)
    real(kind=jprw) :: VDnodes(dstruct%nb_nodes,2)
    real(kind=jprw) :: Lengths(dstruct%nb_nodes)
    real(kind=jprw) :: Qtmp(dstruct%nb_levels,dstruct%nb_nodes)
    real(kind=jprw) :: volD(dstruct%nb_levels,dstruct%nb_nodes)

    real(kind=jprw), pointer :: vol(:,:), S(:,:), pole_bc(:), coords(:,:)

    jlev = 1

    vol     => scalar_field_3d("dual_volumes",dstruct)
    S       => vector_field_2d("dual_normals",dstruct)
    pole_bc => scalar_field_2d("pole_bc",dstruct)
    coords  => vector_field_2d("coordinates",dstruct)

    volD(:,:) = vol(:,:)*D(:,:)
    Lengths(:) = 0.
    VDnodes(:,:) = 0.
    do jedge=1,dstruct%nb_edges
      Sx = S(jedge,XX)
      Sy = S(jedge,YY)
      Ssqr =  Sx*Sx + Sy*Sy 
      V(jedge,XX) = Sx/Ssqr * VDS(jedge)
      V(jedge,YY) = Sy/Ssqr * VDS(jedge)

      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      x1 = coords(ip1,XX)
      x2 = coords(ip2,XX)
      y1 = coords(ip1,YY)
      y2 = coords(ip2,YY)
      length = 1._jprw/sqrt( (x2-x1)**2 + (y2-y1)**2 )
      Lengths(ip1) = Lengths(ip1) + length
      Lengths(ip2) = Lengths(ip2) + length
      VDnodes(ip1,XX) = VDnodes(ip1,XX) + V(jedge,XX) * length
      VDnodes(ip2,XX) = VDnodes(ip2,XX) + V(jedge,XX) * length
      VDnodes(ip1,YY) = VDnodes(ip1,YY) + V(jedge,YY) * length
      VDnodes(ip2,YY) = VDnodes(ip2,YY) + V(jedge,YY) * length
    end do

    call synchronise(VDnodes,dstruct) ! Qmax and Qmin could be synced here
    call synchronise(Lengths,dstruct) ! Qmax and Qmin could be synced here

    do jnode=1,dstruct%nb_nodes
      VDnodes(jnode,:) = VDnodes(jnode,:) / Lengths(jnode)
    end do

    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      V(jedge,:) = 0.5_jprw * (VDnodes(ip1,:) + VDnodes(ip2,:) )
    end do

    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      V(iedge,YY) = 0.
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
      !aun(jedge) = V(jedge,XX)*S(jedge,XX) + V(jedge,YY)*S(jedge,YY)
      aun(jedge) = VDS(jedge)

      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      apos = max(0._jprw,aun(jedge))
      aneg = min(0._jprw,aun(jedge))
      fluxv(jedge) = Q(jlev,ip1)*apos + Q(jlev,ip2)*aneg
    enddo
    !$OMP END PARALLEL DO

     !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,adv,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      adv = 0.0
      if(dstruct%nb_neighbours(jnode) > 1) then
        do jedge = 1,dstruct%nb_neighbours(jnode)
          iedge = dstruct%my_edges(jedge,jnode)
          adv = adv + dstruct%sign(jedge,jnode)*fluxv(iedge)
        enddo
      endif
     ! Update the unknowns in vertices
     Q(jlev,jnode) = Q(jlev,jnode) - adv/max(eps,volD(jlev,jnode)) * dt
     Qtmp(jlev,jnode) = Q(jlev,jnode)
     Q(jlev,jnode) = Q(jlev,jnode)*DR(jlev,jnode)
    enddo
    !$OMP END PARALLEL DO

    call synchronise(Qtmp,dstruct) ! Qmax and Qmin could be synced here
    call synchronise(Q,dstruct) ! Qmax and Qmin could be synced here





    ! 2. Other passes (making the spatial discretisation higher-order)
    ! ----------------------------------------------------------------
    
    do jpass=2,order

      ! Compute derivatives for mpdata
      call compute_gradient(Q, gradQ, Q_is_vector_component, dstruct)

      call synchronise(gradQ,dstruct)

      ! Compute antidiffusive normal velocity in faces

      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,volume_of_two_cells,dQdx,dQdy,Vx,Vy)
      do jedge = 1,dstruct%nb_edges
        ip1 = dstruct%edges(jedge,1)
        ip2 = dstruct%edges(jedge,2)

        ! evaluate gradient and velocity at edge by combining 2 neighbouring dual cells
        volume_of_two_cells = max(eps, volD(jlev,ip1) + volD(jlev,ip2) )
        dQdx = (gradQ(jlev,ip1,XX)+gradQ(jlev,ip2,XX)) / volume_of_two_cells
        dQdy = (gradQ(jlev,ip1,YY)+gradQ(jlev,ip2,YY)) / volume_of_two_cells
        Vx = V(jedge,XX)
        Vy = V(jedge,YY)
        ! variable sign option with asymptotic analysis, (mpdata gauge)
        aun(jedge) = abs(aun(jedge))*(Q(jlev,ip2)-Q(jlev,ip1))*0.5_jprw &
          &          -0.5_jprw*dt*aun(jedge)*(Vx*dQdx+Vy*dQdy)  ! = VDS*dQds
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
            adv = adv + dstruct%sign(jedge,jnode)*aun(iedge)
          enddo
        endif
        ! Update the unknowns in vertices
        Q(jlev,jnode) = Qtmp(jlev,jnode) - adv/max(eps,volD(jlev,jnode)) * dt
        Q(jlev,jnode) = Q(jlev,jnode) * DR(jlev,jnode)
      enddo
      !$OMP END PARALLEL DO

      call synchronise(Q,dstruct)

    end do ! other passes

  contains
    
    subroutine compute_Qmax_and_Qmin( )
      real(kind=jprw) :: Q1, Q2    
      !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge,Q1,Q2)
      do jnode=1,dstruct%nb_nodes
        Q1 = Q(jlev,jnode)
        Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q1 )
        Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q1 )
        do jedge = 1,dstruct%nb_neighbours(jnode)
          Q2 = Q(jlev,dstruct%neighbours(jedge,jnode))
          iedge = dstruct%my_edges(jedge,jnode)
          if (.not. Q_is_vector_component) then
            Qmax(jlev,jnode) = max( Qmax(jlev,jnode), Q2 )
            Qmin(jlev,jnode) = min( Qmin(jlev,jnode), Q2 )
          else
            Qmax(jlev,jnode) = max( Qmax(jlev,jnode), pole_bc(iedge)*Q2 )
            Qmin(jlev,jnode) = min( Qmin(jlev,jnode), pole_bc(iedge)*Q2 )
          end if
        end do
      end do
      !$OMP END PARALLEL DO

      call synchronise(Qmin,dstruct)
      call synchronise(Qmax,dstruct)

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
          apos = max(0._jprw,aun(iedge))
          aneg = min(0._jprw,aun(iedge))
          asignp = max(0._jprw,dstruct%sign(jedge,jnode))
          asignn = min(0._jprw,dstruct%sign(jedge,jnode))
          rhin  = rhin  - asignp*aneg - asignn*apos
          rhout = rhout + asignp*apos + asignn*aneg
        end do
        cp(jnode) = ( Qmax(jlev,jnode)-Q(jlev,jnode) )*volD(jlev,jnode)/( rhin * dt + eps )
        cn(jnode) = ( Q(jlev,jnode)-Qmin(jlev,jnode) )*volD(jlev,jnode)/( rhout* dt + eps )
      end do
      !$OMP END PARALLEL DO

      call synchronise(cp,dstruct)
      call synchronise(cn,dstruct)

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

  end subroutine mpdata_gauge_Q
  
  
  
  
  

  subroutine compute_gradient(Q,gradQ,Q_is_vector_component,dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), intent(in)    :: Q(:,:)
    real(kind=jprw), intent(inout) :: gradQ(:,:,:)
    logical, intent(in) :: Q_is_vector_component
    real(kind=jprw), pointer :: S(:,:)
    real(kind=jprw) :: Sx,Sy,avgQ
    integer :: jedge,iedge,ip1,ip2,jnode,jlev
    real(kind=jprw) :: avgQSx(dstruct%nb_edges)
    real(kind=jprw) :: avgQSy(dstruct%nb_edges)

    jlev = 1

    S   => vector_field_2d("dual_normals",dstruct)

    ! derivatives 

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2,Sx,Sy)
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Sx  = S(jedge,XX)
      Sy  = S(jedge,YY)
      avgQSx(jedge) = Sx*( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw
      avgQSy(jedge) = Sy*( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO SCHEDULE(GUIDED,256) PRIVATE(jnode,jedge,iedge)
    do jnode=1,dstruct%nb_nodes
      gradQ(jlev,jnode,XX) = 0.
      gradQ(jlev,jnode,YY) = 0.
      do jedge = 1,dstruct%nb_neighbours(jnode)
        iedge = dstruct%my_edges(jedge,jnode)
        gradQ(jlev,jnode,XX) = gradQ(jlev,jnode,XX)+dstruct%sign(jedge,jnode)*avgQSx(iedge)
        gradQ(jlev,jnode,YY) = gradQ(jlev,jnode,YY)+dstruct%sign(jedge,jnode)*avgQSy(iedge)
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
        Sy    = S(iedge,YY)
        avgQ  = ( Q(jlev,ip1) + Q(jlev,ip2) )*0.5_jprw

        ! correct for wrong Y-derivatives in previous loop,
        gradQ(jlev,ip2,YY) = gradQ(jlev,ip2,YY) + 2._jprw*Sy*avgQ 

      end do
    end if
  end subroutine compute_gradient

end module mpdata_module


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
  use mpdata_module, &
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

  subroutine setup_isentropic(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw) :: y, cos_y, sin_y
    real(kind=jprw), pointer :: coords(:,:), vol(:,:), hx(:), hy(:), dhxdy_over_G(:,:), pole_bc(:), G(:,:)
    integer :: jnode, jedge, iedge, jlev

    jlev = 1

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
      y = coords(jnode,YY)
      cos_y = cos(y)
      sin_y = sin(y)
      hx(jnode) = radius*cos_y
      hy(jnode) = radius
      G(jlev,jnode) = hx(jnode)*hy(jnode)
      vol(jlev,jnode) = vol(jlev,jnode)*G(jlev,jnode)
      dhxdy_over_G(jlev,jnode) = - sin_y/(radius*max(eps,cos_y))
    enddo

    pole_bc(:) = 1.
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      pole_bc(iedge) = -1.
    end do

    call create_field_in_nodes_2d("coriolis",1,dstruct)
    call create_field_in_nodes_3d("depth",1,dstruct)
    call create_field_in_nodes_3d("velocity",2,dstruct)
    call create_field_in_nodes_3d("velocity_forcing",2,dstruct)
    call create_field_in_nodes_3d("depth_backup",1,dstruct)
    call create_field_in_nodes_3d("velocity_backup",2,dstruct)

    call create_field_in_edges_3d("advective_velocity",2,dstruct)
    call create_field_in_nodes_3d("depth_ratio",1,dstruct) ! old/new
    
    call create_field_in_nodes_2d("topography",1,dstruct)
    call create_field_in_nodes_3d("height",1,dstruct)

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
    real(kind=jprw), pointer :: D(:,:), cor(:), H0(:), H(:,:)
    real(kind=jprw), pointer :: U(:,:,:), coords(:,:)
    integer :: jnode, jlev
    real(kind=jprw) :: x,y
    real(kind=jprw), parameter :: USCAL = 20.
    real(kind=jprw), parameter :: H00 = grav * 8e3
    real(kind=jprw), parameter :: pvel = USCAL/radius
    real(kind=jprw), parameter :: beta = 0.! pi/4._jprw

    jlev = 1

    coords => vector_field_2d("coordinates",dstruct)
    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    H => scalar_field_3d("height",dstruct)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,x,y)
    do jnode=1,dstruct%nb_nodes
      x=coords(jnode,XX)
      y=coords(jnode,YY)
      cor(jnode)   = f0 *( -cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta) )
      D(jlev,jnode)     = (H00-radius**2*(f0+pvel)*0.5*pvel*(-cos(x)*cos(y)*sin(beta)+sin(y)*cos(beta))**2)
      D(jlev,jnode)     = max(0._jprw, D(jlev,jnode)/grav - H0(jnode))
      U(jlev,jnode,XX)  =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y)
      U(jlev,jnode,YY)  = -pvel*sin(x)*sin(beta)*radius
      H(jlev,jnode) = H0(jnode) + D(jlev,jnode)
    end do
    !$OMP END PARALLEL DO

  end subroutine set_state_zonal_flow


  subroutine set_topography(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:), pointer :: H0
    real(kind=jprw), dimension(:,:), pointer :: coords
    real(kind=jprw) :: amp = 2000. ! amplitude of hill
    real(kind=jprw) :: rad = 2.*pi/18. ! radius of hill
    real(kind=jprw) :: xcent = 3.*pi/2.  ! centre of hill
    real(kind=jprw) :: ycent = pi/6.*1.
    real(kind=jprw) :: gamm = 1. ! slope of hill
    real(kind=jprw) :: dist, xlon, ylat
    integer :: jnode

    H0 => scalar_field_2d("topography",dstruct)
    coords => vector_field_2d("coordinates",dstruct)
    
    do jnode=1,dstruct%nb_nodes
      xlon = coords(jnode,XX)
      ylat = coords(jnode,YY)

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
    real(kind=jprw), pointer :: D(:,:), D0(:,:)
    real(kind=jprw), pointer :: Q(:,:,:), Q0(:,:,:)
    integer :: jnode, jlev
    
    jlev = 1

    D  => scalar_field_3d("depth",dstruct)
    D0 => scalar_field_3d("depth_backup",dstruct)
    Q  => vector_field_3d("velocity",dstruct)
    Q0 => vector_field_3d("velocity_backup",dstruct)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      D0(jlev,jnode)   = D(jlev,jnode)
      Q0(jlev,jnode,:) = Q(jlev,jnode,:)
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

    jlev = 1

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
        Rx = R(jlev,jnode,XX)
        Ry = R(jlev,jnode,YY)
        Vnodes(jlev,jnode,XX)=(U(jlev,jnode,XX)+0.5*dt*Rx) 
        Vnodes(jlev,jnode,YY)=(U(jlev,jnode,YY)+0.5*dt*Ry)
      end do
      !$OMP END PARALLEL DO
      call compute_gradient( Vnodes(:,:,XX), grad_Vx, .True., dstruct )
      call compute_gradient( Vnodes(:,:,YY), grad_Vy, .True., dstruct )

      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Vx,Vy,Rx,Ry,dVxdx,dVxdy,dVydx,dVydy)
      do jnode=1,dstruct%nb_nodes
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
      enddo
      !$OMP END PARALLEL DO
      call synchronise( Vnodes, dstruct )
    else if( option .eq. "extrapolate") then
      !dir$ ivdep
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy,U0x,U0y)
      do jnode=1,dstruct%nb_nodes
        Ux    = U(jlev,jnode,XX)
        Uy    = U(jlev,jnode,YY)
        U0x   = U0(jlev,jnode,XX)
        U0y   = U0(jlev,jnode,YY)
        Vnodes(jlev,jnode,XX) = ( 1.5_jprw*Ux - 0.5_jprw*U0x ) * hy(jnode)
        Vnodes(jlev,jnode,YY) = ( 1.5_jprw*Uy - 0.5_jprw*U0y ) * hx(jnode)
      end do
      !$OMP END PARALLEL DO
    end if

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jedge,ip1,ip2)
    do jedge=1,dstruct%nb_edges
      ip1 = dstruct%edges(jedge,1)
      ip2 = dstruct%edges(jedge,2)
      Vedges(jlev,jedge,XX) = (Vnodes(jlev,ip1,XX)+Vnodes(jlev,ip2,XX))*0.5_jprw
      Vedges(jlev,jedge,YY) = (Vnodes(jlev,ip1,YY)+Vnodes(jlev,ip2,YY))*0.5_jprw
    enddo
    !$OMP END PARALLEL DO
    
    ! Since the pole point lies outside the lon-lat domain, Vedges is wrongly calculated
    ! y_pole .ne. 0.5(y1+y2)
    do jedge=1,dstruct%nb_pole_edges
      iedge = dstruct%pole_edges(jedge)
      Vedges(jlev,iedge,YY) = 0.
    enddo

  end subroutine compute_advective_velocities

  subroutine compute_forcing(dstruct)
    type(DataStructure_type), intent(inout) :: dstruct
    integer :: jnode, jlev
    real(kind=jprw) :: Ux, Uy
    real(kind=jprw), dimension(:),   pointer :: H0, hx, hy, cor
    real(kind=jprw), dimension(:,:), pointer :: H, D, vol, dhxdy_over_G, coords
    real(kind=jprw), dimension(:,:,:), pointer :: U, R
    real(kind=jprw) :: grad_H(dstruct%nb_levels,dstruct%nb_nodes, 2)

    jlev = 1

    coords => vector_field_2d("coordinates",dstruct)
    vol => scalar_field_3d("dual_volumes",dstruct)
    hx => scalar_field_2d("hx",dstruct)
    hy => scalar_field_2d("hy",dstruct)
    dhxdy_over_G => scalar_field_3d("dhxdy_over_G",dstruct)
    cor => scalar_field_2d("coriolis",dstruct)

    H  => scalar_field_3d("height",dstruct)
    H0 => scalar_field_2d("topography",dstruct)
    D  => scalar_field_3d("depth",dstruct)
    U  => vector_field_3d("velocity",dstruct)
    R  => vector_field_3d("velocity_forcing",dstruct)

    do jlev=1,dstruct%nb_levels
      H(jlev,:) = H0(:) + D(jlev,:)
    end do
    jlev = 1

    call compute_gradient( H, grad_H, .False., dstruct )

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy)
    do jnode=1,dstruct%nb_nodes
      Ux    = U(jlev,jnode,XX)
      Uy    = U(jlev,jnode,YY)
      R(jlev,jnode,XX) = -grav*grad_H(jlev,jnode,XX)*hy(jnode)/vol(jlev,jnode) &
        &           + cor(jnode)*Uy - dhxdy_over_G(jlev,jnode)*Ux*Uy
      R(jlev,jnode,YY) = -grav*grad_H(jlev,jnode,YY)*hx(jnode)/vol(jlev,jnode) &
        &           - cor(jnode)*Ux + dhxdy_over_G(jlev,jnode)*Ux*Ux
    end do
    !$OMP END PARALLEL DO

    call synchronise(R,dstruct)
  end subroutine compute_forcing



  subroutine add_forcing_to_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), pointer :: U(:,:,:) , R(:,:,:), D(:,:)
    integer :: jnode, jlev

    jlev = 1

    D => scalar_field_3d("depth",dstruct)
    U => vector_field_3d("velocity",dstruct)
    R => vector_field_3d("velocity_forcing",dstruct)
    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode)
    do jnode=1,dstruct%nb_nodes
      U(jlev,jnode,XX) = ( U(jlev,jnode,XX) + 0.5_jprw*dt*R(jlev,jnode,XX) )
      U(jlev,jnode,YY) = ( U(jlev,jnode,YY) + 0.5_jprw*dt*R(jlev,jnode,YY) )  
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

    jlev = 1 
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
    jlev = 1

    call compute_gradient( H, grad_H, .False., dstruct )

    !dir$ ivdep
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jnode,Ux,Uy,Rx_exp,Ry_exp,Ux_adv,Uy_adv,m,Rx,Ry)
    do jnode=1,dstruct%nb_nodes
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
    !$OMP END PARALLEL DO
    call synchronise(U,dstruct)
    call synchronise(R,dstruct)

  end subroutine implicit_solve



  subroutine advect_solution(dt,dstruct)
    real(kind=jprw), intent(in) :: dt
    type(DataStructure_type), intent(inout) :: dstruct
    real(kind=jprw), dimension(:,:),   pointer :: D, D0, DR
    real(kind=jprw), dimension(:,:,:), pointer :: U, V
    real(kind=jprw) :: VDS(dstruct%nb_edges)
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
