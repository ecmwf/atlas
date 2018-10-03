! (C) Copyright 2013 ECMWF.
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

module atlas_nabla_program_module
use atlas_module
use, intrinsic :: iso_c_binding, only: c_double, c_int
implicit none

  character(len=1024) :: msg

  type(atlas_Grid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_functionspace_NodeColumns) :: node_columns
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_fvm_Method) :: fvm
  type(atlas_Nabla) :: nabla
  type(atlas_Field) :: varfield
  type(atlas_Field) :: gradfield
  type(atlas_Config) :: config

  real(c_double), pointer :: var(:,:)
  real(c_double), pointer :: grad(:,:,:)

  integer :: nlev
  integer :: niter
  integer :: nouter
  integer :: startcount


  integer, parameter :: JPRB = c_double
  integer, parameter :: JPIM = c_int
  real(JPRB), parameter :: RPI = 2.0_JPRB*asin(1.0_JPRB)
  real(JPRB) :: RA

  type :: timer_type
  private
    integer*8 :: clck_counts_start, clck_counts_stop, clck_rate
    integer*8 :: counted = 0
    logical :: paused = .True.
  contains
    procedure, public :: start   => timer_start
    procedure, public :: pause   => timer_pause
    procedure, public :: resume  => timer_resume
    procedure, public :: elapsed => timer_elapsed
  end type timer_type

contains

    function timer_elapsed(self) result(time)
        use, intrinsic :: iso_c_binding, only : c_double
        class(timer_type), intent(inout) :: self
        real(c_double) :: time
        if (.not. self%paused) then
            call system_clock ( self%clck_counts_stop, self%clck_rate )
            time = (self%counted + self%clck_counts_stop - self%clck_counts_start)/real(self%clck_rate)
        else if (self%counted .ge. 0) then
            time = self%counted/real(self%clck_rate)
        else
            time = 0.
        end if
      end function timer_elapsed

    subroutine timer_start(self)
        class(timer_type), intent(inout) :: self
        call system_clock ( self%clck_counts_start, self%clck_rate )
        self%paused = .False.
        self%counted = 0
    end subroutine timer_start

    subroutine timer_pause(self)
        class(timer_type), intent(inout) :: self
        call system_clock ( self%clck_counts_stop, self%clck_rate )
        self%counted = self%counted + self%clck_counts_stop - self%clck_counts_start
        self%paused = .True.
    end subroutine timer_pause

    subroutine timer_resume(self)
        class(timer_type), intent(inout) :: self
        call system_clock ( self%clck_counts_start, self%clck_rate )
        self%paused = .False.
    end subroutine timer_resume

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

TYPE(ATLAS_FIELD) :: LONLAT,DUAL_VOLUMES,DUAL_NORMALS,NODE2EDGE_SIGN
type(atlas_mesh_Nodes) :: NODES
type(atlas_mesh_Edges) :: EDGES
TYPE(ATLAS_CONNECTIVITY) :: EDGE2NODE, NODE2EDGE
REAL(KIND=JPRB), POINTER :: ZLONLAT(:,:),ZDUAL_VOLUMES(:),ZDUAL_NORMALS(:,:),&
 & ZNODE2EDGE_SIGN(:,:)
INTEGER(KIND=ATLAS_KIND_IDX),POINTER :: IEDGE2NODE(:,:)
INTEGER(KIND=ATLAS_KIND_IDX),POINTER :: INODE2EDGE(:)
INTEGER(KIND=JPIM) :: INODE2EDGE_SIZE
INTEGER(KIND=JPIM) :: JNODE,JEDGE,JLEV,INEDGES,IP1,IP2,IEDGE,INODES
REAL(KIND=JPRB) :: ZAVG,ZSIGN,ZMETRIC_X,ZMETRIC_Y,ZSCALE
REAL(KIND=JPRB), ALLOCATABLE :: ZAVG_S(:,:,:)
!     ------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('FV_GRADIENT',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>nlev,&
 & node_columns=>fvm )

!write(0,*) 'enter fv_gradient'
!write(0,*) 'shape pvar ',shape(pvar)
NODES = MESH%NODES()
EDGES = MESH%EDGES()
LONLAT = NODES%LONLAT()
EDGE2NODE = EDGES%NODE_CONNECTIVITY()
NODE2EDGE = NODES%EDGE_CONNECTIVITY()
DUAL_VOLUMES=NODES%FIELD('dual_volumes')
DUAL_NORMALS=EDGES%FIELD('dual_normals')
NODE2EDGE_SIGN=NODES%FIELD('node2edge_sign')

CALL LONLAT%DATA(ZLONLAT)
CALL EDGE2NODE%DATA(IEDGE2NODE)
CALL DUAL_VOLUMES%DATA(ZDUAL_VOLUMES)
CALL DUAL_NORMALS%DATA(ZDUAL_NORMALS)
CALL NODE2EDGE_SIGN%DATA(ZNODE2EDGE_SIGN)

INEDGES = SIZE(ZDUAL_NORMALS)/2
ZSCALE = RPI/180.0_JPRB * RPI/180.0_JPRB * RA

ALLOCATE(ZAVG_S(2,NFLEVG,INEDGES))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JEDGE,IP1,IP2,JLEV,ZAVG)
DO JEDGE=1,INEDGES
  IP1 = IEDGE2NODE(1,JEDGE)
  IP2 = IEDGE2NODE(2,JEDGE)
  DO JLEV=1,NFLEVG
    ZAVG = (PVAR(JLEV,IP1)+PVAR(JLEV,IP2)) * 0.5_JPRB
    ZAVG_S(1,JLEV,JEDGE) = ZDUAL_NORMALS(1,JEDGE)*RPI/180.0_JPRB * ZAVG
    ZAVG_S(2,JLEV,JEDGE) = ZDUAL_NORMALS(2,JEDGE)*RPI/180.0_JPRB * ZAVG
  ENDDO
ENDDO
!$OMP END PARALLEL DO

INODES = SIZE(ZLONLAT)/2
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JNODE,JLEV,JEDGE,IEDGE,ZSIGN,ZMETRIC_X,ZMETRIC_Y,INODE2EDGE,INODE2EDGE_SIZE)
DO JNODE=1,INODES
  DO JLEV=1,NFLEVG
    PGRAD(1,JLEV,JNODE) = 0.0
    PGRAD(2,JLEV,JNODE) = 0.0
  ENDDO
  CALL NODE2EDGE%ROW(JNODE,INODE2EDGE,INODE2EDGE_SIZE)
  DO JEDGE=1,INODE2EDGE_SIZE
    IEDGE = INODE2EDGE(JEDGE)
    ZSIGN = ZNODE2EDGE_SIGN(JEDGE,JNODE)
    DO JLEV=1,NFLEVG
      PGRAD(1,JLEV,JNODE) =  PGRAD(1,JLEV,JNODE)+ZSIGN*ZAVG_S(1,JLEV,IEDGE)
      PGRAD(2,JLEV,JNODE) =  PGRAD(2,JLEV,JNODE)+ZSIGN*ZAVG_S(2,JLEV,IEDGE)
    ENDDO
  ENDDO
  ZMETRIC_Y = 1._JPRB/(ZDUAL_VOLUMES(JNODE)*ZSCALE)
  ZMETRIC_X = ZMETRIC_Y/COS(ZLONLAT(2,JNODE)*RPI/180.0_JPRB)
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
  use fckit_resource_module
  character(len=:),allocatable :: grid_uid
  type(atlas_mesh_Nodes) :: mesh_nodes

  call atlas_library%initialise()

  call fckit_resource("--grid","N24",grid_uid)
  call fckit_resource("--levels",137,nlev)
  call fckit_resource("--iterations",100,niter)
  call fckit_resource("--startcount",5,startcount)
  call fckit_resource("--outer",1,nouter)

  config = atlas_Config()
  call config%set("radius",1.0)
  if( .not.config%get("radius",RA) ) RA = 1.0

  ! Setup
  grid = atlas_Grid(grid_uid)
  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid) ! second optional argument for atlas_GridDistrubution
  fvm  = atlas_fvm_Method(mesh,config)
  node_columns = fvm%node_columns()
  nabla = atlas_Nabla(fvm)

  ! Create a variable field and a gradient field
  varfield   = node_columns%create_field(name="var", kind=atlas_real(c_double), levels=nlev)
  gradfield  = node_columns%create_field(name="grad", kind=atlas_real(c_double),levels=nlev, variables=2)

  ! Access to data
  call varfield%data(var)
  call gradfield%data(grad)
  var(:,:) = 0.

  mesh_nodes = mesh%nodes()
  write(msg,*) "Mesh has locally ",mesh_nodes%size(), " nodes"
  call atlas_log%info(msg)
  call mesh_nodes%final()

end subroutine

! -----------------------------------------------------------------------------

subroutine finalize()
  ! Cleanup
  call config%final()
  call varfield%final()
  call gradfield%final()
  call nabla%final()
  call fvm%final()
  call node_columns%final()
  call nodes%final()
  call mesh%final()
  call grid%final()
  call meshgenerator%final()
  call atlas_library%finalise()
end subroutine

! -----------------------------------------------------------------------------

subroutine run()
use fckit_mpi_module
type(timer_type) :: timer;
type(fckit_mpi_comm) :: mpi
integer :: jiter, jouter
real(c_double) :: timing_cpp, timing_f90, timing
real(c_double) :: min_timing_cpp, min_timing_f90
mpi = fckit_mpi_comm()
call node_columns%halo_exchange(varfield)
call mpi%barrier()
timing_cpp = 1.e10
timing_f90 = 1.e10
min_timing_cpp = 1.e10
min_timing_f90 = 1.e10

do jouter=1,nouter
! Compute the gradient
timing = 1.e10
do jiter = 1,niter
    call timer%start()
    call nabla%gradient(varfield,gradfield)
    timing = min(timing,timer%elapsed())
enddo
timing_cpp = timing
write(msg,*) "timing_cpp = ", timing_cpp
call atlas_log%info(msg)
if( nouter == 1 .or. jouter < nouter ) then
  min_timing_cpp = min(timing_cpp,min_timing_cpp)
endif

! Compute the gradient with Fortran routine above
timing = 1.e10
do jiter = 1,niter
  call timer%start()
  call FV_GRADIENT(var,grad)
  timing = min(timing,timer%elapsed())
enddo
timing_f90 = timing
write(msg,*) "timing_f90 = ", timing_f90
call atlas_log%info(msg)
if( nouter == 1 .or. jouter < nouter ) then
  min_timing_f90 = min(timing_f90,min_timing_f90)
endif
write(msg,*) "|timing_f90-timing_cpp| / timing_f90 = ", &
  & abs(timing_f90-timing_cpp)/timing_f90 *100 , "%"
call atlas_log%info(msg)

enddo

call atlas_log%info("==================")
write(msg,*) "min_timing_cpp = ", min_timing_cpp
call atlas_log%info(msg)
write(msg,*) "min_timing_f90 = ", min_timing_f90
call atlas_log%info(msg)
write(msg,*) "|min_timing_f90-min_timing_cpp| / min_timing_f90 = ", &
  & abs(min_timing_f90-min_timing_cpp)/min_timing_f90 *100 , "%"
call atlas_log%info(msg)
end subroutine

end module atlas_nabla_program_module

! -----------------------------------------------------------------------------

program atlas_nabla_program
use atlas_nabla_program_module
  call init()
  call run()
  call finalize()
end program
