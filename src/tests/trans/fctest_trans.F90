! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_trans_fixture
use atlas_module
use fckit_module
use iso_c_binding
implicit none
character(len=1024) :: msg
type(fckit_mpi_comm) :: mpi
end module fctest_atlas_trans_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_trans,fctest_atlas_trans_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call fckit_main%init()
  mpi = fckit_mpi_comm()
  if( mpi%rank() /= 0 ) then
    call fckit_log%reset()
  endif
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  write(msg,*) "TESTSUITE_FINALIZE"; call fckit_log%info(msg)
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_trans_backend )
  type(atlas_Trans) :: trans
  FCTEST_CHECK( trans%has_backend("local") )
  FCTEST_CHECK( trans%has_backend("ifs") )
  if( trans%has_backend("ifs") ) then
    FCTEST_CHECK_EQUAL( trans%backend(), "ifs" )
  else
    FCTEST_CHECK_EQUAL( trans%backend(), "local" )
  endif  
END_TEST

! -----------------------------------------------------------------------------

TEST( test_trans )
  type(atlas_StructuredGrid) :: grid
  type(atlas_StructuredGrid) :: trans_grid
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Mesh) :: mesh
  type(atlas_Trans) :: trans
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_functionspace_NodeColumns) :: nodes_fs
  type(atlas_functionspace_Spectral) :: spectral_fs, trans_spectral
  type(atlas_Field)         :: scalarfield1, scalarfield2
  type(atlas_Field)         :: windfield
  type(atlas_Field)         :: vorfield,divfield
  type(atlas_Field)         :: spectralfield1, spectralfield2
  type(atlas_FieldSet)      :: scalarfields
  type(atlas_FieldSet)      :: spectralfields
  real(c_double), pointer :: scal1(:,:), scal2(:), spec1(:,:), spec2(:), wind(:,:,:), vor(:,:), div(:,:)
  real(c_double), allocatable :: check(:)
  integer :: nlev, truncation, jn, in, jlev
  integer, pointer :: nvalue(:)
  type(atlas_Field) :: glb_vorfield

  real(c_double) :: tol

  tol = 1.e-8
  nlev=10
  truncation = 21

  grid = atlas_StructuredGrid("O24")

  FCTEST_CHECK_EQUAL( grid%owners(), 1 )

  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid)

  call meshgenerator%final()

  FCTEST_CHECK_EQUAL( mesh%owners(), 1 )
  FCTEST_CHECK_EQUAL( grid%owners(), 2 ) ! mesh tracks grid

  trans = atlas_Trans(grid,truncation)
  FCTEST_CHECK_EQUAL( grid%owners(), 3 ) ! trans tracks grid

  ! FCTEST_CHECK_EQUAL( trans%nb_gridpoints_global(), int(grid%size()) )

  trans_grid = trans%grid()
  FCTEST_CHECK_EQUAL( trans_grid%owners(), 4 )

  FCTEST_CHECK( .not. trans%is_null() )
  FCTEST_CHECK_EQUAL( trans%truncation(), truncation )
  
  trans_spectral = trans%spectral()
  FCTEST_CHECK_EQUAL( trans_spectral%nb_spectral_coefficients_global(), (truncation+1)*(truncation+2) )

  nodes = mesh%nodes()
  nodes_fs = atlas_functionspace_NodeColumns(mesh,0)

  write(msg,*) "nodes_fs%owners()",nodes_fs%owners(); call fckit_log%info(msg)

  scalarfield1 = nodes_fs%create_field(name="scalar1",kind=atlas_real(c_double),levels=nlev)
  write(msg,*) "nodes_fs%owners()",nodes_fs%owners(); call fckit_log%info(msg)

  scalarfield2 = nodes_fs%create_field(name="scalar2",kind=atlas_real(c_double))
  write(msg,*) "nodes_fs%owners()",nodes_fs%owners(); call fckit_log%info(msg)

  spectral_fs = atlas_functionspace_Spectral(trans)
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)

  FCTEST_CHECK_EQUAL( spectral_fs%nb_spectral_coefficients_global(), (truncation+1)*(truncation+2) )

  spectralfield1 = spectral_fs%create_field(name="spectral1",kind=atlas_real(c_double),levels=nlev)
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)

  spectralfield2 = spectral_fs%create_field(name="spectral2",kind=atlas_real(c_double))
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)

  call scalarfield1%data(scal1)
  call scalarfield2%data(scal2)
  call spectralfield1%data(spec1)
  call spectralfield2%data(spec2)

  write(msg,*) "shape = ", spectralfield2%shape(); call fckit_log%info(msg)

  ! All waves to zero except wave 1 to 3
  spec1(1:nlev,:) = 0
  spec1(1:nlev,1) = 3
  ! All waves to zero except wave 1 to 4
  spec2(:) = 0
  spec2(1) = 4

  call trans%invtrans(spectralfield1,scalarfield1)
  call trans%dirtrans(scalarfield1,spectralfield1)

  allocate( check(nlev) )
  check(:) = 3
  FCTEST_CHECK_CLOSE( spec1(:,1), check, tol )
  deallocate( check )

  scalarfields = atlas_FieldSet("scalarfields")
  call scalarfields%add(scalarfield1)
  call scalarfields%add(scalarfield2)

  spectralfields = atlas_FieldSet("spectralfields")
  call spectralfields%add(spectralfield1)
  call spectralfields%add(spectralfield2)

  call trans%invtrans(spectralfields,scalarfields)
  call trans%dirtrans(scalarfields,spectralfields)

  allocate( check(nlev) )
  check(:) = 3
  FCTEST_CHECK_CLOSE( spec1(:,1), check, tol )
  check(:) = 0
  FCTEST_CHECK_CLOSE( spec1(:,2), check, tol )
  FCTEST_CHECK_CLOSE( spec1(:,3), check, tol )
  FCTEST_CHECK_CLOSE( spec1(:,4), check, tol )
  FCTEST_CHECK_CLOSE( spec1(:,5), check, tol )
  deallocate( check )

  FCTEST_CHECK_CLOSE( spec2(1), 4._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(2), 0._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(3), 0._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(4), 0._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(5), 0._c_double, tol )

  windfield = nodes_fs%create_field(name="wind",kind=atlas_real(c_double),levels=nlev,variables=3)
  call windfield%data(wind)
  wind(:,:,:) = 0._c_double
  write(msg,*) "nodes_fs%owners()",nodes_fs%owners(); call fckit_log%info(msg)

  vorfield = spectral_fs%create_field(name="vorticity",kind=atlas_real(c_double),levels=nlev)
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners()

  call vorfield%data(vor)

  divfield =  spectral_fs%create_field(name="divergence",kind=atlas_real(c_double),levels=nlev)
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)

  call divfield%data(div)

  call trans%dirtrans_wind2vordiv(windfield,vorfield,divfield)

  call trans%invtrans_vordiv2wind(vorfield,divfield,windfield)

  glb_vorfield = spectral_fs%create_field(name="vorticity",kind=atlas_real(c_double),levels=nlev,global=.true.)
  call spectral_fs%gather(vorfield,glb_vorfield)
  call spectral_fs%scatter(glb_vorfield,vorfield)

  write(msg,*) "cleaning up"; call fckit_log%info(msg)
  write(msg,*) "nodes_fs%owners()",nodes_fs%owners(); call fckit_log%info(msg)
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)

  call scalarfield1%final()
  call scalarfield2%final()
  call spectralfield1%final()
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)
  call spectralfield2%final()
  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)
  call windfield%final()
  call vorfield%final()
  call divfield%final()
  call glb_vorfield%final()

  write(msg,*) "nodes_fs%owners()",nodes_fs%owners(); call fckit_log%info(msg)
  call nodes_fs%final()

  write(msg,*) "spectral_fs%owners()",spectral_fs%owners(); call fckit_log%info(msg)
  call spectral_fs%final()

  call scalarfields%final()
  call spectralfields%final()
  call mesh%final()
  call trans%final()
  call grid%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_trans_nomesh )
  type(atlas_StructuredGrid) :: grid
  type(atlas_Trans) :: trans
  type(atlas_Partitioner) :: partitioner
  type(atlas_functionspace_StructuredColumns) :: gridpoints_fs
  type(atlas_functionspace_Spectral) :: spectral_fs
  type(atlas_Field)         :: scalarfield1, scalarfield2
  type(atlas_Field)         :: spectralfield1, spectralfield2
  type(atlas_FieldSet)      :: scalarfields
  type(atlas_FieldSet)      :: spectralfields
  real(c_double), pointer :: scal1(:,:), scal2(:), spec1(:,:), spec2(:)
  real(c_double), allocatable :: check(:)
  integer :: nlev, truncation, jn, in, jlev
  integer, pointer :: nvalue(:)
  real(c_double) :: tol

  tol = 1.e-8
  nlev=10
  truncation = 21

  grid = atlas_StructuredGrid("O24")
  trans = atlas_Trans(grid,truncation)

  partitioner = atlas_Partitioner(type="ifs")
  gridpoints_fs = atlas_functionspace_StructuredColumns(grid,partitioner)
  scalarfield1 = gridpoints_fs%create_field(name="scalar1",kind=atlas_real(c_double),levels=nlev)
  scalarfield2 = gridpoints_fs%create_field(name="scalar2",kind=atlas_real(c_double))

  spectral_fs = atlas_functionspace_Spectral(trans)
  spectralfield1 = spectral_fs%create_field(name="spectral1",kind=atlas_real(c_double),levels=nlev)
  spectralfield2 = spectral_fs%create_field(name="spectral2",kind=atlas_real(c_double))

  call scalarfield1%data(scal1)
  call scalarfield2%data(scal2)
  call spectralfield1%data(spec1)
  call spectralfield2%data(spec2)

  ! All waves to zero except wave 1 to 3
  spec1(1:nlev,:) = 0
  spec1(1:nlev,1) = 3
  ! All waves to zero except wave 1 to 4
  spec2(:) = 0
  spec2(1) = 4

  call atlas_log%debug("invtrans")
  call trans%invtrans(spectralfield1,scalarfield1)
  call atlas_log%debug("dirtrans")
  call trans%dirtrans(scalarfield1,spectralfield1)

  allocate( check(nlev) )
  check(:) = 3
  FCTEST_CHECK_CLOSE( spec1(:,1), check, tol )
  deallocate( check )

  scalarfields = atlas_FieldSet("scalarfields")
  call scalarfields%add(scalarfield1)
  call scalarfields%add(scalarfield2)

  spectralfields = atlas_FieldSet("spectralfields")
  call spectralfields%add(spectralfield1)
  call spectralfields%add(spectralfield2)

  call trans%invtrans(spectralfields,scalarfields)
  call trans%dirtrans(scalarfields,spectralfields)

  allocate( check(nlev) )
  check(:) = 3
  FCTEST_CHECK_CLOSE( spec1(:,1), check, tol )
  check(:) = 0
  FCTEST_CHECK_CLOSE( spec1(:,2), check, tol )
  FCTEST_CHECK_CLOSE( spec1(:,3), check, tol )
  FCTEST_CHECK_CLOSE( spec1(:,4), check, tol )
  FCTEST_CHECK_CLOSE( spec1(:,5), check, tol )
  deallocate( check )

  FCTEST_CHECK_CLOSE( spec2(1), 4._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(2), 0._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(3), 0._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(4), 0._c_double, tol )
  FCTEST_CHECK_CLOSE( spec2(5), 0._c_double, tol )

  write(msg,*) "cleaning up"; call fckit_log%info(msg)

  call scalarfield1%final()
  call scalarfield2%final()
  call spectralfield1%final()
  call spectralfield2%final()
  call spectral_fs%final()
  call gridpoints_fs%final()
  call scalarfields%final()
  call spectralfields%final()
  call trans%final()
  call grid%final()
END_TEST

TEST( test_transdwarf )
type(atlas_StructuredGrid) :: grid
type(atlas_Partitioner) :: partitioner
type(atlas_Trans) :: trans
type(atlas_functionspace_StructuredColumns) :: gridpoints
type(atlas_functionspace_Spectral) :: spectral
type(atlas_Field) :: fieldg, field
type(atlas_FieldSet) :: gpfields, spfields
integer :: jfld, nfld
character(len=10) :: fieldname
real(c_double) :: norm
real(c_double), pointer :: gvar(:)

grid = atlas_StructuredGrid("O24")
trans = atlas_Trans(grid,23)
partitioner = atlas_Partitioner("ifs")
gridpoints = atlas_functionspace_StructuredColumns(grid,partitioner)
spectral = atlas_functionspace_Spectral(trans)

gpfields = atlas_FieldSet("gridpoint")
spfields = atlas_FieldSet("spectral")

nfld=10
do jfld=1,nfld
  fieldg = gridpoints%create_field(atlas_real(c_double),global=.true.)
  field  = gridpoints%create_field(atlas_real(c_double))

  ! Read global field data
  ! ...
  FCTEST_CHECK_EQUAL( fieldg%rank(), 1 )
  call fieldg%data(gvar)
  if( mpi%rank() == 0 ) then
    gvar(:) = 0.
  endif

  call gridpoints%scatter(fieldg,field)

  call gpfields%add( field )
  call spfields%add( spectral%create_field(atlas_real(c_double)) )

  FCTEST_CHECK_EQUAL( field%owners(), 2 )
enddo

call trans%dirtrans(gpfields,spfields)
call trans%invtrans(spfields,gpfields)

do jfld=1,spfields%size()
  field = spfields%field(jfld)
  write(msg,*) "spectral field ",field%name(); call atlas_log%info(msg)
  call spectral%norm(field,norm)
  write(msg,*) "norm = ",norm; call atlas_log%info(msg)
enddo

call field%final()
call fieldg%final()
call gpfields%final()
call spfields%final()
call gridpoints%final()
call spectral%final()
call trans%final()
call grid%final()
END_TEST

TEST( test_spectral_only )
type(atlas_functionspace_Spectral) :: spectral
type(atlas_Field) :: field, fieldg
integer :: jfld, nfld, T
character(len=10) :: fieldname
real(c_double) :: norm
integer, parameter :: RE=0
integer, parameter :: IM=1
integer :: jm, jn, n, m, jc
real(c_double), pointer :: sp(:,:)
integer, pointer :: nmyms(:), nvalue(:), nasm0(:)


spectral = atlas_functionspace_Spectral(truncation=21,levels=10)
T = spectral%truncation()
nfld = spectral%levels()

field = spectral%create_field(kind=atlas_real(c_double),levels=0)
FCTEST_CHECK_EQUAL(field%rank(), 1)

field = spectral%create_field(atlas_real(c_double))
FCTEST_CHECK_EQUAL(field%rank(), 2)
FCTEST_CHECK_EQUAL(field%shape(1), nfld)

! Initialise distributed field, plus checks
    call field%data(sp)
    nmyms  => spectral%nmyms()
    nvalue => spectral%nvalue()
    nasm0  => spectral%nasm0()
    jc = 1
    do jm=1,spectral%nump()
      m = nmyms(jm)
      FCTEST_CHECK_EQUAL( nasm0(m), jc )
      do jn=0,T-m
        n = m + jn
        FCTEST_CHECK( jc <= spectral%nb_spectral_coefficients() - 1 )
        FCTEST_CHECK_EQUAL( nvalue(jc), n )
        do jfld=1,nfld
          sp(jfld,jc+RE) = sp_value(m,n,RE)
          sp(jfld,jc+IM) = sp_value(m,n,IM)
        enddo
        jc = jc+2
      enddo
    enddo

fieldg = spectral%create_field(atlas_real(c_double),global=.true.)
call spectral%gather(field,fieldg)

! Check global field (gathered matches distributed field)
    call fieldg%data(sp)
    T = spectral%truncation()
    if( mpi%rank() == 0 ) then
        write(0,*) spectral%nb_spectral_coefficients_global()
        jc = 1
        do m=0,T
          do n=m,T
            FCTEST_CHECK( jc <= spectral%nb_spectral_coefficients_global() - 1 )
            do jfld=1,nfld
              FCTEST_CHECK_EQUAL( sp(jfld,jc+RE), sp_value(m,n,RE) )
              FCTEST_CHECK_EQUAL( sp(jfld,jc+IM), sp_value(m,n,IM) )
            enddo
            jc = jc+2
          enddo
        enddo
    endif

call field%final()
call spectral%final()

contains
function sp_value(m,n,complex_component)
  real(c_double) :: sp_value
  integer :: m, n, complex_component
  if( complex_component == RE ) then
    sp_value = m*T+n
  else
    if( m == 0 ) then
      sp_value = 0
    else
      sp_value = -m*T+n
    endif
  endif
end function

END_TEST


END_TESTSUITE

