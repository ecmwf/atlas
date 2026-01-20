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

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_fvm_nabla_Fxt
use atlas_module
use, intrinsic :: iso_c_binding
implicit none

  type(atlas_StructuredGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_fvm_Method) :: fvm
  type(atlas_Nabla) :: nabla
  type(atlas_functionspace_NodeColumns) :: node_columns
  type(atlas_Field) :: varfield
  type(atlas_Field) :: gradfield

  real(c_double), pointer :: var(:,:)
  real(c_double), pointer :: grad(:,:,:)

  integer, parameter :: nlev = 137

  type(atlas_Config) :: config

  integer, parameter :: JPRB = c_double
  integer, parameter :: JPIM = c_int
  real(JPRB), parameter :: RPI = 2.0_JPRB*asin(1.0_JPRB)
  real(JPRB) :: RA



  type :: timer_type
  private
    integer(8) :: clck_counts_start, clck_counts_stop, clck_rate
    integer(8) :: counted = 0
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
            time = real(self%counted + self%clck_counts_stop - self%clck_counts_start,c_double) &
                & /real(self%clck_rate,c_double)
        else if (self%counted .ge. 0) then
            time = real(self%counted,c_double)/real(self%clck_rate,c_double)
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


    subroutine rotated_flow_magnitude( fvm, field, beta, radius )
        type(atlas_fvm_method), intent(in) :: fvm
        type(atlas_field), intent(inout) :: field
        real(c_double), intent(in) :: beta
        real(c_double), intent(in) :: radius

        !real(c_double) :: radius
        real(c_double) :: USCAL
        real(c_double) :: pvel
        real(c_double) :: deg2rad

        type(atlas_functionspace_NodeColumns) :: fs_nodecolumns
        type(atlas_mesh_Nodes) :: nodes
        type(atlas_field) :: lonlat
        real(c_double), pointer :: lonlat_deg(:,:)
        real(c_double), pointer :: var(:,:)
        integer :: nnodes, jnode, nlev, jlev
        real(c_double) :: x, y, Ux, Uy
        integer :: LON=1
        integer :: LAT=2

        !radius = fvm%radius()
        USCAL = 20.
        pvel = USCAL/radius
        deg2rad = RPI/180.

        fs_nodecolumns = fvm%node_columns()
        mesh = fs_nodecolumns%mesh()
        nodes = fs_nodecolumns%nodes()
        lonlat = nodes%lonlat()
        call lonlat%data(lonlat_deg)
        call field%data(var)
        nnodes = nodes%size()
        nlev = field%levels()

        do jnode=1,nnodes
           x = lonlat_deg(LON,jnode) * deg2rad
           y = lonlat_deg(LAT,jnode) * deg2rad
           Ux =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y)
           Uy = -pvel*sin(x)*sin(beta)*radius

           do jlev=1,nlev
             var(jlev,jnode) = sqrt(Ux*Ux+Uy*Uy)
           enddo

        enddo



    end subroutine


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
REAL(KIND=JPRB) :: ZAVG,ZSIGN,ZMETRIC_X,ZMETRIC_Y
REAL(KIND=JPRB), ALLOCATABLE :: ZAVG_S(:,:,:)
INTEGER(KIND=JPIM) :: NFLEVG
REAL(KIND=JPRB) :: SCALE, DEG2RAD
DEG2RAD = RPI/180.
SCALE = DEG2RAD*DEG2RAD*RA

!     ------------------------------------------------------------------

!IF (LHOOK) CALL DR_HOOK('FV_GRADIENT',0,ZHOOK_HANDLE)

NFLEVG=nlev

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
ALLOCATE(ZAVG_S(2,NFLEVG,INEDGES))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JEDGE,IP1,IP2,JLEV,ZAVG)
DO JEDGE=1,INEDGES
  IP1 = IEDGE2NODE(1,JEDGE)
  IP2 = IEDGE2NODE(2,JEDGE)
  DO JLEV=1,NFLEVG
    ZAVG = (PVAR(JLEV,IP1)+PVAR(JLEV,IP2)) * 0.5_JPRB
    ZAVG_S(1,JLEV,JEDGE) = ZDUAL_NORMALS(1,JEDGE)*DEG2RAD*ZAVG
    ZAVG_S(2,JLEV,JEDGE) = ZDUAL_NORMALS(2,JEDGE)*DEG2RAD*ZAVG
  ENDDO
ENDDO
!$OMP END PARALLEL DO

INODES = SIZE(ZLONLAT)/2
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(JNODE,JLEV,JEDGE,IEDGE,ZSIGN,ZMETRIC_X,ZMETRIC_Y,INODE2EDGE,INODE2EDGE_SIZE)
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
  ZMETRIC_Y = 1./(ZDUAL_VOLUMES(JNODE)*SCALE)
  ZMETRIC_X = ZMETRIC_Y/COS(ZLONLAT(2,JNODE)*DEG2RAD)
  DO JLEV=1,NFLEVG
    PGRAD(1,JLEV,JNODE) = PGRAD(1,JLEV,JNODE)*ZMETRIC_X
    PGRAD(2,JLEV,JNODE) = PGRAD(2,JLEV,JNODE)*ZMETRIC_Y
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!IF (LHOOK) CALL DR_HOOK('FV_GRADIENT',1,ZHOOK_HANDLE)
END SUBROUTINE FV_GRADIENT




end module fctest_atlas_fvm_nabla_Fxt

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_fvm_nabla,fctest_atlas_fvm_nabla_Fxt)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()

  RA = 10.

  config = atlas_Config()
  call config%set("radius",RA)
  call config%set("levels",nlev)

  ! Setup
  grid = atlas_StructuredGrid("N24")
  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid) ! second optional argument for atlas_GridDistrubution
  fvm  = atlas_fvm_Method(mesh,config)
  node_columns = fvm%node_columns()
  nabla = atlas_Nabla(fvm)

  ! Create a variable field and a gradient field
  varfield  = node_columns%create_field(name="var",kind=atlas_real(c_double),levels=nlev)
  gradfield = node_columns%create_field(name="grad",kind=atlas_real(c_double),levels=nlev,variables=2)

  ! Access to data
  call varfield%data(var)
  call gradfield%data(grad)
  var(:,:) = 0.

  call rotated_flow_magnitude(fvm,varfield,beta=0.5*RPI*0.75,radius=RA)


END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
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
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_fvm )
type(atlas_StructuredGrid) :: grid
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_Mesh) :: mesh
type(atlas_fvm_Method) :: fvm

grid = atlas_StructuredGrid("N24")
meshgenerator = atlas_MeshGenerator()
mesh = meshgenerator%generate(grid)
fvm  = atlas_fvm_Method(mesh)

call fvm%final()
call mesh%final()
call grid%final()

END_TEST

! -----------------------------------------------------------------------------

TEST( test_nabla )
type(timer_type) :: timer;
integer :: jiter, niter
real(c_double), allocatable :: norm_native(:)
real(c_double), allocatable :: norm_fortran(:)
real(c_double) :: checked_value_X = 1.7893197319163668E-016
real(c_double) :: checked_value_Y = -1.0547670904632068E-006

type(atlas_Output) :: gmsh

call node_columns%halo_exchange(varfield)

niter = 5

! Compute the gradient
call timer%start()
do jiter = 1,niter
call nabla%gradient(varfield,gradfield)
enddo

write(0,*) "time elapsed: ", timer%elapsed()
call node_columns%halo_exchange(gradfield)
call node_columns%mean(gradfield,norm_native)
write(0,*) "mean : ", norm_native

FCTEST_CHECK_CLOSE( norm_native(1), checked_value_X, 1.e-12_c_double )
FCTEST_CHECK_CLOSE( norm_native(2), checked_value_Y, 1.e-12_c_double )

write(0,*) ""

! Compute the gradient with Fortran routine above
call timer%start()
do jiter = 1,niter
CALL FV_GRADIENT(var,grad)
enddo

write(0,*) "time elapsed: ", timer%elapsed()
call node_columns%halo_exchange(gradfield)
call node_columns%mean(gradfield,norm_fortran)
write(0,*) "mean : ", norm_fortran

FCTEST_CHECK_CLOSE( norm_native(1), checked_value_X, 1.e-12_c_double )
FCTEST_CHECK_CLOSE( norm_native(2), checked_value_Y, 1.e-12_c_double )

gmsh = atlas_output_Gmsh("out_atlas_fctest_fvm_nabla.msh", levels=[10])
call gmsh%write( mesh )
call gmsh%write( varfield )
call gmsh%write( gradfield )
END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

