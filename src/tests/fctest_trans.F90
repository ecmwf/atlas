! (C) Copyright 1996-2016 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_trans_fixture
use atlas_module
use iso_c_binding
implicit none
end module fctest_atlas_trans_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_trans,fctest_atlas_trans_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_mpi_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_trans )
  type(atlas_ReducedGrid) :: grid
  type(atlas_Mesh) :: mesh
  type(atlas_Trans) :: trans
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_functionspace_Nodes) :: nodes_fs
  type(atlas_functionspace_Spectral) :: spectral_fs
  type(atlas_Field)         :: scalarfield1, scalarfield2
  type(atlas_Field)         :: windfield
  type(atlas_Field)         :: vorfield,divfield
  type(atlas_Field)         :: spectralfield1, spectralfield2
  type(atlas_FieldSet)      :: scalarfields
  type(atlas_FieldSet)      :: spectralfields
  real(c_double), pointer :: scal1(:,:), scal2(:), spec1(:,:), spec2(:), wind(:,:,:), vor(:,:), div(:,:)
  real(c_double), allocatable :: check(:)
  integer :: nlev, nsmax, jn, in, jlev
  integer, pointer :: nvalue(:)
  type(atlas_Field) :: glb_vorfield

  real(c_double) :: tol

  tol = 1.e-8
  nlev=10
  nsmax = 21

  grid = atlas_ReducedGrid("O24")

  FCTEST_CHECK_EQUAL( grid%owners(), 1 )

  mesh = atlas_generate_mesh(grid)

  FCTEST_CHECK_EQUAL( mesh%owners(), 1 )

  trans = atlas_Trans(grid,nsmax)

  FCTEST_CHECK( .not. trans%is_null() )
  !call atlas_write_gmsh(mesh,"testf3.msh")
  !FCTEST_CHECK_EQUAL( trans%owners(), 1 )
  !FCTEST_CHECK_EQUAL( trans%owners(), 1 )
!  FCTEST_CHECK_EQUAL( trans%owners(), 1 )
!  FCTEST_CHECK_EQUAL( trans%owners(), 1 )
  FCTEST_CHECK_EQUAL( trans%nproc(), 1 )
  FCTEST_CHECK_EQUAL( trans%myproc(proc0=1), 1 )
  FCTEST_CHECK_EQUAL( trans%ndgl(), grid%nlat() )
  FCTEST_CHECK_EQUAL( trans%ngptot(), grid%npts() )
  FCTEST_CHECK_EQUAL( trans%ngptotg(), grid%npts() )
  FCTEST_CHECK_EQUAL( trans%nsmax(), nsmax )

  nodes = mesh%nodes()
  nodes_fs = atlas_functionspace_Nodes(mesh,0)

  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  scalarfield1 = nodes_fs%create_field("scalar1",atlas_real(c_double),nlev)
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  scalarfield2 = nodes_fs%create_field("scalar2",atlas_real(c_double))
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  spectral_fs = atlas_functionspace_Spectral(trans)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  spectralfield1 = spectral_fs%create_field("spectral1",nlev)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  spectralfield2 = spectral_fs%create_field("spectral2")
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call scalarfield1%data(scal1)
  call scalarfield2%data(scal2)
  call spectralfield1%data(spec1)
  call spectralfield2%data(spec2)

  write(0,*) "shape = ", spectralfield2%shape()

  ! All waves to zero except wave 1 to 3
  spec1(1:nlev,:) = 0
  spec1(1:nlev,1) = 3
  ! All waves to zero except wave 1 to 4
  spec2(:) = 0
  spec2(1) = 4

  call trans%invtrans(spectral_fs,spectralfield1,nodes_fs,scalarfield1)
  call trans%dirtrans(nodes_fs,scalarfield1,spectral_fs,spectralfield1)

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

  call trans%invtrans(spectral_fs,spectralfields,nodes_fs,scalarfields)
  call trans%dirtrans(nodes_fs,scalarfields,spectral_fs,spectralfields)

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

  windfield = nodes_fs%create_field("wind",atlas_real(c_double),nlev,(/3/))
  call windfield%data(wind)
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  vorfield = spectral_fs%create_field("vorticity",nlev)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call vorfield%data(vor)

  divfield =  spectral_fs%create_field("divergence",nlev)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call divfield%data(div)

  call trans%dirtrans_wind2vordiv(nodes_fs,windfield,spectral_fs,vorfield,divfield)

  nvalue => trans%nvalue()
  FCTEST_CHECK_EQUAL( size(nvalue), trans%nspec2() )

  do jlev=1,nlev
    do jn=1,trans%nspec2()
      in = nvalue(jn)
      if( in > 5 ) then
        vor(jlev,jn) = 0
      endif
    enddo
  enddo

  call trans%invtrans_vordiv2wind(spectral_fs,vorfield,divfield,nodes_fs,windfield)

  glb_vorfield = spectral_fs%create_global_field("vorticity",nlev)
  call spectral_fs%gather(vorfield,glb_vorfield)
  call spectral_fs%scatter(glb_vorfield,vorfield)

  write(0,*) "cleaning up"
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call scalarfield1%final()
  call scalarfield2%final()
  call spectralfield1%final()
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()
  call spectralfield2%final()
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()
  call windfield%final()
  call vorfield%final()
  call divfield%final()
  call glb_vorfield%final()

  write(0,*) "nodes_fs%owners()",nodes_fs%owners()
  call nodes_fs%final()

  write(0,*) "spectral_fs%owners()",spectral_fs%owners()
  call spectral_fs%final()

  call scalarfields%final()
  call spectralfields%final()
  call mesh%final()
  call trans%final()
  call grid%final()
END_TEST

! -----------------------------------------------------------------------------

TEST( test_trans_nomesh )
  type(atlas_ReducedGrid) :: grid
  type(atlas_Trans) :: trans
  type(atlas_functionspace_ReducedGridPoint) :: gridpoints_fs
  type(atlas_functionspace_Spectral) :: spectral_fs
  type(atlas_Field)         :: scalarfield1, scalarfield2
  type(atlas_Field)         :: spectralfield1, spectralfield2
  type(atlas_FieldSet)      :: scalarfields
  type(atlas_FieldSet)      :: spectralfields
  real(c_double), pointer :: scal1(:,:), scal2(:), spec1(:,:), spec2(:)
  real(c_double), allocatable :: check(:)
  integer :: nlev, nsmax, jn, in, jlev
  integer, pointer :: nvalue(:)
  real(c_double) :: tol

  tol = 1.e-8
  nlev=10
  nsmax = 21

  grid = atlas_ReducedGrid("O24")
  trans = atlas_Trans(grid,nsmax)

  gridpoints_fs = atlas_functionspace_ReducedGridPoint(grid)
  scalarfield1 = gridpoints_fs%create_field("scalar1",nlev)
  scalarfield2 = gridpoints_fs%create_field("scalar2")

  spectral_fs = atlas_functionspace_Spectral(trans)
  spectralfield1 = spectral_fs%create_field("spectral1",nlev)
  spectralfield2 = spectral_fs%create_field("spectral2")

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

  write(0,*) "cleaning up"

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
type(atlas_ReducedGrid) :: grid
type(atlas_Trans) :: trans
type(atlas_functionspace_ReducedGridPoint) :: gridpoints
type(atlas_functionspace_Spectral) :: spectral
type(atlas_Field) :: fieldg, field
type(atlas_FieldSet) :: gpfields, spfields
integer :: jfld, nfld
character(len=10) :: fieldname

grid = atlas_ReducedGrid("O24")
trans = atlas_Trans(grid,23)
gridpoints = atlas_functionspace_ReducedGridPoint(grid)
spectral = atlas_functionspace_Spectral(trans)

gpfields = atlas_FieldSet("gridpoint")
spfields = atlas_FieldSet("spectral")

nfld=10
do jfld=1,nfld
  write(fieldname,'(I0)') jfld

  fieldg = gridpoints%create_global_field(fieldname)
  field  = gridpoints%create_field(fieldname)

  ! Read global field data
  ! ...

  call gridpoints%scatter(fieldg,field)

  call gpfields%add( field )
  call spfields%add( spectral%create_field(fieldname) )

  FCTEST_CHECK_EQUAL( field%owners(), 2 )
enddo

call trans%dirtrans(gpfields,spfields)
call trans%invtrans(spfields,gpfields)

do jfld=1,spfields%size()
  field = spfields%field(jfld)
  write(atlas_log%msg,*) "spectral field ",field%name(); call atlas_log%info()
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


END_TESTSUITE

