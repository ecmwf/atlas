! (C) Copyright 1996-2015 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! gradient computations using Edge-Based-Finite-Volume method
! A comparison between C++ built-in and a fortran version is done here
! @author Willem Deconinck

! -----------------------------------------------------------------------------

module atlas_numerics_nabla_module
use atlas_module
use iso_c_binding, only: c_double, c_int
implicit none

  type(atlas_ReducedGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_functionspace_EdgeBasedFiniteVolume) :: fvm
  type(atlas_Nabla) :: nabla
  type(atlas_Field) :: varfield
  type(atlas_Field) :: gradfield

  real(c_double), pointer :: var(:,:)
  real(c_double), pointer :: grad(:,:,:)

  integer :: nlev

  type(atlas_Config) :: config

  integer, parameter :: JPRB = c_double
  integer, parameter :: JPIM = c_int
  real(JPRB), parameter :: RPI = 2.0_JPRB*asin(1.0_JPRB)
  real(JPRB) :: RA
  type(atlas_FunctionSpace) :: edges

  type :: Timer_type
  private
    integer*8 :: clck_counts_start, clck_counts_stop, clck_rate
    integer*8 :: counted = 0
    logical :: paused = .True.
  contains
    procedure, public :: start   => Timer_start
    procedure, public :: pause   => Timer_pause
    procedure, public :: resume  => Timer_resume
    procedure, public :: elapsed => Timer_elapsed
  end type Timer_type

contains

    function Timer_elapsed(self) result(time)
        use iso_c_binding, only : c_double
        class(Timer_type), intent(inout) :: self
        real(c_double) :: time
        if (.not. self%paused) then
            call system_clock ( self%clck_counts_stop, self%clck_rate )
            time = (self%counted + self%clck_counts_stop - self%clck_counts_start)/real(self%clck_rate)
        else if (self%counted .ge. 0) then
            time = self%counted/real(self%clck_rate)
        else
            time = 0.
        end if
      end function Timer_elapsed

    subroutine Timer_start(self)
        class(Timer_type), intent(inout) :: self
        call system_clock ( self%clck_counts_start, self%clck_rate )
        self%paused = .False.
        self%counted = 0
    end subroutine Timer_start

    subroutine Timer_pause(self)
        class(Timer_type), intent(inout) :: self
        call system_clock ( self%clck_counts_stop, self%clck_rate )
        self%counted = self%counted + self%clck_counts_stop - self%clck_counts_start
        self%paused = .True.
    end subroutine Timer_pause

    subroutine Timer_resume(self)
        class(Timer_type), intent(inout) :: self
        call system_clock ( self%clck_counts_start, self%clck_rate )
        self%paused = .False.
    end subroutine Timer_resume

SUBROUTINE FV_GRADIENT(PVAR,PGRAD)

!**** *GP_DERIVATIVES* - COMPUTES GRID-POINT DERIVATIVES using FINITE VOLUME

!     PURPOSE.

!**   INTERFACE.
!     ----------
!        *CALL* *GP_DERIVATIVES(..)*

!        EXPLICIT ARGUMENTS
!        --------------------


!        IMPLICIT ARGUMENTS :
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.     NONE.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      Mats Hamrud AND Willem Deconinck  *ECMWF*
!      ORIGINAL : 88-02-04

!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

!USE PARKIND1 , ONLY : JPIM     ,JPRB
!USE ATLAS_MODULE
!USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK
!USE YOMGEM   , ONLY : TGEM
!USE YOMATLAS , ONLY : YRATLAS
!USE YOMDIM   , ONLY : YRDIM
!USE YOMDIMV  , ONLY : YRDIMV
!USE YOMCST   , ONLY : RPI, RA

IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN)   :: PVAR(:,:)
REAL(KIND=JPRB),INTENT(OUT)  :: PGRAD(:,:,:)

TYPE(ATLAS_FIELD) :: LONLAT,EDGE2NODE,DUAL_VOLUMES,DUAL_NORMALS,&
 NODE2EDGE,NODE2EDGE_SIZE,NODE2EDGE_SIGN
TYPE(ATLAS_MESH_NODES) :: NODES

REAL(KIND=JPRB), POINTER :: ZLONLAT(:,:),ZDUAL_VOLUMES(:),ZDUAL_NORMALS(:,:),&
 & ZNODE2EDGE_SIGN(:,:)
INTEGER(KIND=JPIM),POINTER :: IEDGE2NODE(:,:),INODE2EDGE(:,:),INODE2EDGE_SIZE(:)
INTEGER(KIND=JPIM) :: JNODE,JEDGE,JLEV,INEDGES,IP1,IP2,IEDGE,INODES
REAL(KIND=JPRB) :: ZAVG,ZSIGN,ZMETRIC_X,ZMETRIC_Y
REAL(KIND=JPRB), ALLOCATABLE :: ZAVG_S(:,:,:)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('FV_GRADIENT',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>nlev,&
 & NODES_FS=>fvm )

!write(0,*) 'enter fv_gradient'
!write(0,*) 'shape pvar ',shape(pvar)
NODES = MESH%NODES()
LONLAT = NODES%LONLAT()
EDGE2NODE=EDGES%FIELD('nodes')
DUAL_VOLUMES=NODES%FIELD('dual_volumes')
DUAL_NORMALS=EDGES%FIELD('dual_normals')
NODE2EDGE=NODES%FIELD('to_edge')
NODE2EDGE_SIZE=NODES%FIELD('to_edge_size')
NODE2EDGE_SIGN=NODES%FIELD('node2edge_sign')

CALL LONLAT%DATA(ZLONLAT)
CALL EDGE2NODE%DATA(IEDGE2NODE)
CALL DUAL_VOLUMES%DATA(ZDUAL_VOLUMES)
CALL DUAL_NORMALS%DATA(ZDUAL_NORMALS)
CALL NODE2EDGE%DATA(INODE2EDGE)
CALL NODE2EDGE_SIZE%DATA(INODE2EDGE_SIZE)
CALL NODE2EDGE_SIGN%DATA(ZNODE2EDGE_SIGN)

INEDGES = SIZE(ZDUAL_NORMALS)/2
ALLOCATE(ZAVG_S(2,NFLEVG,INEDGES))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JEDGE,IP1,IP2,JLEV,ZAVG)
DO JEDGE=1,INEDGES
  IP1 = IEDGE2NODE(1,JEDGE)
  IP2 = IEDGE2NODE(2,JEDGE)
  DO JLEV=1,NFLEVG
    ZAVG = (PVAR(JLEV,IP1)+PVAR(JLEV,IP2)) * 0.5_JPRB
    ZAVG_S(1,JLEV,JEDGE) = ZDUAL_NORMALS(1,JEDGE)*ZAVG
    ZAVG_S(2,JLEV,JEDGE) = ZDUAL_NORMALS(2,JEDGE)*ZAVG
  ENDDO
ENDDO
!$OMP END PARALLEL DO

INODES = SIZE(ZLONLAT)/2
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JNODE,JLEV,JEDGE,IEDGE,ZSIGN,ZMETRIC_X,ZMETRIC_Y)
DO JNODE=1,INODES
  DO JLEV=1,NFLEVG
    PGRAD(1,JLEV,JNODE) = 0.0
    PGRAD(2,JLEV,JNODE) = 0.0
  ENDDO
  DO JEDGE=1,INODE2EDGE_SIZE(JNODE)
    IEDGE = INODE2EDGE(JEDGE,JNODE)
    ZSIGN = ZNODE2EDGE_SIGN(JEDGE,JNODE)
    DO JLEV=1,NFLEVG
      PGRAD(1,JLEV,JNODE) =  PGRAD(1,JLEV,JNODE)+ZSIGN*ZAVG_S(1,JLEV,IEDGE)
      PGRAD(2,JLEV,JNODE) =  PGRAD(2,JLEV,JNODE)+ZSIGN*ZAVG_S(2,JLEV,IEDGE)
    ENDDO
  ENDDO
  ZMETRIC_X = RA/ZDUAL_VOLUMES(JNODE)
  ZMETRIC_Y = ZMETRIC_X*COS(ZLONLAT(2,JNODE)*RPI/180.0)
  DO JLEV=1,NFLEVG
    PGRAD(1,JLEV,JNODE) = PGRAD(1,JLEV,JNODE)*ZMETRIC_X
    PGRAD(2,JLEV,JNODE) = PGRAD(2,JLEV,JNODE)*ZMETRIC_Y
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END ASSOCIATE

!IF (LHOOK) CALL DR_HOOK('FV_GRADIENT',1,ZHOOK_HANDLE)
END SUBROUTINE FV_GRADIENT



subroutine init()
  character(len=1024) :: grid_uid

  call atlas_init()

  call atlas_resource("--grid","N24",grid_uid)
  call atlas_resource("--levels",137,nlev)

  config = atlas_Config()
  call config%set("radius",1.0)

  ! Setup
  grid = atlas_ReducedGrid(grid_uid)
  meshgenerator = atlas_ReducedGridMeshGenerator()
  mesh = meshgenerator%generate(grid) ! second optional argument for atlas_GridDistrubution
  fvm  = atlas_functionspace_EdgeBasedFiniteVolume(mesh,config)
  nabla = atlas_Nabla(fvm)

  ! Create a variable field and a gradient field
  varfield = fvm%create_field("var",atlas_real(c_double),nlev)
  gradfield  = fvm%create_field("grad",atlas_real(c_double),nlev,[2])

  ! Access to data
  call varfield%data(var)
  call gradfield%data(grad)
  var(:,:) = 0.

  if( .not.config%get("radius",RA) ) RA = 1.0
  edges = mesh%function_space("edges")

end subroutine

! -----------------------------------------------------------------------------

subroutine finalize()
  ! Cleanup
  call config%final()
  call varfield%final()
  call gradfield%final()
  call nabla%final()
  call fvm%final()
  call nodes%final()
  call mesh%final()
  call grid%final()
  call meshgenerator%final()
  call atlas_finalize()
end subroutine

! -----------------------------------------------------------------------------

subroutine run()
type(Timer_type) :: timer
integer :: jiter, niter
real(c_double) :: timing_cpp, timing_f90

call fvm%halo_exchange(varfield)

niter = 5

! Compute the gradient
call timer%start()
do jiter = 1,niter
call nabla%gradient(varfield,gradfield)
enddo
timing_cpp = timer%elapsed()
write(0,*) "timing_cpp = ", timing_cpp

! Compute the gradient with Fortran routine above
call timer%start()
do jiter = 1,niter
CALL FV_GRADIENT(var,grad)
enddo
timing_f90 = timer%elapsed()
write(0,*) "timing_f90 = ", timing_f90

write(0,*) "|timing_f90-timing_cpp| / timing_f90 = ", abs(timing_f90-timing_cpp)/timing_f90 *100 , "%"

end subroutine

end module atlas_numerics_nabla_module

! -----------------------------------------------------------------------------

program atlas_numerics_nabla
use atlas_numerics_nabla_module
  call init()
  call run()
  call finalize()
end program
