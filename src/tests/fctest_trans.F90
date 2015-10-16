! (C) Copyright 1996-2015 ECMWF.
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
  type(atlas_Nodes) :: nodes
  type(atlas_NodesFunctionSpace) :: nodes_fs
  type(atlas_SpectralFunctionspace) :: spectral_fs
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
  real(c_double), allocatable :: vorg(:,:)

  real(c_double) :: tol

  tol = 1.e-8
  nlev=10
  nsmax = 21

  grid = atlas_ReducedGrid("oct.N24")

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
  nodes_fs = atlas_NodesFunctionSpace("nodes",mesh,0)

  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  scalarfield1 = nodes_fs%create_field("scalar1",atlas_real(c_double),nlev)
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  scalarfield2 = nodes_fs%create_field("scalar2",atlas_real(c_double))
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  spectral_fs = atlas_SpectralFunctionSpace("spectral",trans)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  spectralfield1 = spectral_fs%create_field("spectral1",nlev)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  spectralfield2 = spectral_fs%create_field("spectral2")
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call scalarfield1%access_data(scal1)
  call scalarfield2%access_data(scal2)
  call spectralfield1%access_data(spec1)
  call spectralfield2%access_data(spec2)

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
  call scalarfields%add_field(scalarfield1)
  call scalarfields%add_field(scalarfield2)

  spectralfields = atlas_FieldSet("spectralfields")
  call spectralfields%add_field(spectralfield1)
  call spectralfields%add_field(spectralfield2)

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
  call windfield%access_data(wind)
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()

  vorfield = spectral_fs%create_field("vorticity",nlev)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call vorfield%access_data(vor)

  divfield =  spectral_fs%create_field("divergence",nlev)
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call divfield%access_data(div)

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

  allocate( vorg( nlev, trans%nspec2g() ) )
  call trans%gathspec(vor,vorg)

  write(0,*) "cleaning up"
  write(0,*) "nodes_fs%owners()",nodes_fs%owners()
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()

  call scalarfield1%finalize()
  call scalarfield2%finalize()
  call spectralfield1%finalize()
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()
  call spectralfield2%finalize()
  write(0,*) "spectral_fs%owners()",spectral_fs%owners()
  call windfield%finalize()
  call vorfield%finalize()
  call divfield%finalize()

  write(0,*) "nodes_fs%owners()",nodes_fs%owners()
  call nodes_fs%finalize()

  write(0,*) "spectral_fs%owners()",spectral_fs%owners()
  call spectral_fs%finalize()

  call scalarfields%finalize()
  call spectralfields%finalize()
  call mesh%finalize()
  call trans%finalize()
  call grid%finalize()
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

