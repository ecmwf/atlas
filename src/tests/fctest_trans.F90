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
use atlas_grids_module
use iso_c_binding
implicit none
end module fctest_atlas_trans_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_trans,fctest_atlas_trans_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
  call atlas_grids_load()
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
  type(atlas_Functionspace) :: nodes
  type(atlas_Functionspace) :: spectral
  type(atlas_Field)         :: scalarfield1, scalarfield2
  type(atlas_Field)         :: windfield
  type(atlas_Field)         :: vorfield,divfield
  type(atlas_Field)         :: spectralfield1, spectralfield2
  type(atlas_FieldSet)      :: scalarfields
  type(atlas_FieldSet)      :: spectralfields
  real(c_double), pointer :: scal1(:,:), scal2(:,:), spec1(:,:), spec2(:,:), wind(:,:), vor(:), div(:)
  real(c_double), allocatable :: check(:)
  integer :: nlev, nsmax, jn, in
  integer, pointer :: nvalue(:)

  real(c_double) :: tol

  tol = 1.e-12
  nlev=100
  nsmax = 21

  grid = atlas_ReducedGrid("oct.N24")
  mesh = atlas_generate_mesh(grid)
  trans = atlas_Trans(grid,nsmax)

  !call atlas_write_gmsh(mesh,"testf3.msh")

  FCTEST_CHECK_EQUAL( trans%nproc(), 1 )
  FCTEST_CHECK_EQUAL( trans%myproc(proc0=1), 1 )
  FCTEST_CHECK_EQUAL( trans%ndgl(), 24*2 )
  FCTEST_CHECK_EQUAL( trans%ngptot(), grid%npts() )
  FCTEST_CHECK_EQUAL( trans%ngptotg(), grid%npts() )
  FCTEST_CHECK_EQUAL( trans%nsmax(), nsmax )

  nodes = mesh%function_space("nodes")
  call nodes%create_field("scalar1",nlev)
  scalarfield1 = nodes%field("scalar1")
  call nodes%create_field("scalar2",1)
  scalarfield2 = nodes%field("scalar2")

  ! note this is in C++ ordering for now (TODO fix)
  call mesh%create_function_space("spectral","spectral",(/trans%nspec2(),ATLAS_FIELD_NB_VARS/))
  spectral = mesh%function_space("spectral")
  call spectral%create_field("spectral1",nlev)
  spectralfield1 = spectral%field("spectral1")
  call spectral%create_field("spectral2",1)
  spectralfield2 = spectral%field("spectral2")

  call scalarfield1%access_data(scal1)
  call scalarfield2%access_data(scal2)
  call spectralfield1%access_data(spec1)
  call spectralfield2%access_data(spec2)

  spec1(:,:) = 0
  spec1(:,1) = 3
  spec2(:,:) = 0
  spec2(:,1) = 4

  call trans%invtrans(spectralfield1,scalarfield1)
  call trans%dirtrans(scalarfield1,spectralfield1)

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

  allocate( check(1) )
  check(:) = 4
  FCTEST_CHECK_CLOSE( spec2(:,1), check, tol )
  check(:) = 0
  FCTEST_CHECK_CLOSE( spec2(:,2), check, tol )
  FCTEST_CHECK_CLOSE( spec2(:,3), check, tol )
  FCTEST_CHECK_CLOSE( spec2(:,4), check, tol )
  FCTEST_CHECK_CLOSE( spec2(:,5), check, tol )
  deallocate( check )


  call nodes%create_field("wind",3)
  windfield = nodes%field("wind")
  call windfield%access_data(wind)

  call spectral%create_field("vorticity",1)
  vorfield = spectral%field("vorticity")
  call vorfield%access_data(vor)

  call spectral%create_field("divergence",1)
  divfield = spectral%field("divergence")
  call vorfield%access_data(div)

  call trans%dirtrans_wind2vordiv(windfield,vorfield,divfield)

  nvalue => trans%nvalue()
  FCTEST_CHECK_EQUAL( size(nvalue), trans%nspec2() )

  do jn=1,trans%nspec2()
    in = nvalue(jn)
    if( in > 20 ) then
      vor(jn) = 0
    endif
  enddo

  call trans%invtrans_vordiv2wind(vorfield,divfield,windfield)

  call atlas_delete(scalarfields)
  call atlas_delete(spectralfields)
  call atlas_delete(mesh)
  call atlas_delete(trans)
  call atlas_delete(grid)
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

