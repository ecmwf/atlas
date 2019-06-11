! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------


TESTSUITE(fctest_atlas_trans_unstr)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  use atlas_module
  type(atlas_Trans) :: trans
  call atlas_library%initialise()
  call trans%set_backend("local")
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use atlas_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_trans )
  use atlas_module
  implicit none
  type(atlas_UnstructuredGrid) :: grid
  type(atlas_Trans) :: trans
  type(atlas_functionspace_Spectral) :: spectral
  type(atlas_Field)         :: gp_scal_field, gp_wind_field
  type(atlas_Field)         :: sp_scal_field, sp_vor_field, sp_div_field
  type(atlas_Config) :: config
  integer :: truncation

  real(c_double) :: xy(2,10)

  xy(:,1) = [107.196,26.0167]
  xy(:,2) = [144.768,28.3721]
  xy(:,3) = [150.891,60.0343]
  xy(:,4) = [164.566,25.5053]
  xy(:,5) = [116.851,14.0295]
  xy(:,6) = [124.84,28.3978]
  xy(:,7) = [157.901,42.042]
  xy(:,8) = [111.41,43.1056]
  xy(:,9) = [134.333,44.6677]
  xy(:,10) = [120.307,59.7167]

  grid = atlas_UnstructuredGrid( xy )

  truncation = 1279

  trans = atlas_Trans(grid,truncation)

  FCTEST_CHECK_EQUAL( trans%truncation(), truncation )

  spectral = trans%spectral()

  sp_scal_field = spectral%create_field(name="spectral_scalar",kind=atlas_real(8))
  sp_vor_field  = spectral%create_field(name="spectral_vorticity",kind=atlas_real(8))
  sp_div_field  = spectral%create_field(name="spectral_divergence",kind=atlas_real(8))
  
  gp_scal_field = atlas_Field(name="gridpoint_scalar",kind=atlas_real(8),shape=[10])
  gp_wind_field = atlas_Field(name="gridpoint_wind",kind=atlas_real(8),shape=[10,2])

  config = atlas_Config()
  call config%set("warning",0) ! turn off warnings for unstructured grids

  call trans%invtrans( sp_scal_field, gp_scal_field, config )
  call trans%invtrans_vordiv2wind( sp_vor_field,  sp_div_field, gp_wind_field, config )

END_TEST


END_TESTSUITE

