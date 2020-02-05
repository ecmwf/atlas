! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Connectivity
! @author Willem Deconinck

#include "fckit/fctest.h"


! -----------------------------------------------------------------------------

TESTSUITE(fctest_atlas_Connectivity)

! -----------------------------------------------------------------------------
TESTSUITE_INIT()
  use fckit_module
  call fckit_main%init()
END_TESTSUITE_INIT

TEST( test_connectivity )
#if 1
  use fckit_module
  use atlas_connectivity_module
  use atlas_kinds_module
  use, intrinsic :: iso_c_binding

  implicit none
  type(atlas_Connectivity) :: connectivity
  integer(ATLAS_KIND_IDX), pointer :: padded(:,:), row(:), data(:,:)
  integer(ATLAS_KIND_IDX), pointer :: cols(:)
  integer(c_int) :: ncols

  call fckit_log%info("test_connectivity starting")

  connectivity = atlas_Connectivity("hybrid")
  FCTEST_CHECK_EQUAL( connectivity%owners(), 1 )

  FCTEST_CHECK_EQUAL(connectivity%name(),"hybrid")

  FCTEST_CHECK_EQUAL(connectivity%rows(),0)
  FCTEST_CHECK_EQUAL(connectivity%missing_value(),0)

  call connectivity%add(2,4, &
    & [ 1, 2, 3, 4,  &
    &   5, 6, 7, 8 ] )

  FCTEST_CHECK_EQUAL(connectivity%mincols(),4)
  FCTEST_CHECK_EQUAL(connectivity%maxcols(),4)
  FCTEST_CHECK_EQUAL(connectivity%rows(),   2)


  call connectivity%data(data,ncols)
  FCTEST_CHECK_EQUAL(ncols,4)
  FCTEST_CHECK_EQUAL(data(1,1), 1)
  FCTEST_CHECK_EQUAL(data(2,1), 2)
  FCTEST_CHECK_EQUAL(data(3,1), 3)
  FCTEST_CHECK_EQUAL(data(4,1), 4)
  FCTEST_CHECK_EQUAL(data(1,2), 5)
  FCTEST_CHECK_EQUAL(data(2,2), 6)
  FCTEST_CHECK_EQUAL(data(3,2), 7)
  FCTEST_CHECK_EQUAL(data(4,2), 8)

  call connectivity%add(2,3, &
    & [ 9,  10, 11,  &
    &   12, 13, 14 ] )

  FCTEST_CHECK_EQUAL(connectivity%mincols(),3)
  FCTEST_CHECK_EQUAL(connectivity%maxcols(),4)
  FCTEST_CHECK_EQUAL(connectivity%rows(),   4)

  call connectivity%add(2,4)
  call connectivity%add(2,3)

  !============= Functional access =============!

  FCTEST_CHECK_EQUAL(connectivity%value(1,1), 1)
  FCTEST_CHECK_EQUAL(connectivity%value(2,1), 2)
  FCTEST_CHECK_EQUAL(connectivity%value(3,1), 3)
  FCTEST_CHECK_EQUAL(connectivity%value(4,1), 4)
  FCTEST_CHECK_EQUAL(connectivity%value(1,2), 5)
  FCTEST_CHECK_EQUAL(connectivity%value(2,2), 6)
  FCTEST_CHECK_EQUAL(connectivity%value(3,2), 7)
  FCTEST_CHECK_EQUAL(connectivity%value(4,2), 8)

  FCTEST_CHECK_EQUAL(connectivity%value(1,3), 9)
  FCTEST_CHECK_EQUAL(connectivity%value(2,3), 10)
  FCTEST_CHECK_EQUAL(connectivity%value(3,3), 11)
  FCTEST_CHECK_EQUAL(connectivity%value(1,4), 12)
  FCTEST_CHECK_EQUAL(connectivity%value(2,4), 13)
  FCTEST_CHECK_EQUAL(connectivity%value(3,4), 14)

  !============= Padded Data pointer access =============!


  ! Todo: Check this aborts,
  !       because data has become irregular
  ! call connectivity%data(data,ncols)


  call connectivity%padded_data(padded,cols)
  FCTEST_CHECK_EQUAL(padded(1,1), 1)
  FCTEST_CHECK_EQUAL(padded(2,1), 2)
  FCTEST_CHECK_EQUAL(padded(3,1), 3)
  FCTEST_CHECK_EQUAL(padded(4,1), 4)
  FCTEST_CHECK_EQUAL(padded(1,2), 5)
  FCTEST_CHECK_EQUAL(padded(2,2), 6)
  FCTEST_CHECK_EQUAL(padded(3,2), 7)
  FCTEST_CHECK_EQUAL(padded(4,2), 8)

  FCTEST_CHECK_EQUAL(padded(1,3), 9)
  FCTEST_CHECK_EQUAL(padded(2,3), 10)
  FCTEST_CHECK_EQUAL(padded(3,3), 11)
  FCTEST_CHECK_EQUAL(padded(1,4), 12)
  FCTEST_CHECK_EQUAL(padded(2,4), 13)
  FCTEST_CHECK_EQUAL(padded(3,4), 14)

  !============= Structural access =============!

  FCTEST_CHECK_EQUAL(connectivity%access%row(1)%col(1), 1)
  FCTEST_CHECK_EQUAL(connectivity%access%row(1)%col(2), 2)
  FCTEST_CHECK_EQUAL(connectivity%access%row(1)%col(3), 3)
  FCTEST_CHECK_EQUAL(connectivity%access%row(1)%col(4), 4)
  FCTEST_CHECK_EQUAL(connectivity%access%row(2)%col(1), 5)
  FCTEST_CHECK_EQUAL(connectivity%access%row(2)%col(2), 6)
  FCTEST_CHECK_EQUAL(connectivity%access%row(2)%col(3), 7)
  FCTEST_CHECK_EQUAL(connectivity%access%row(2)%col(4), 8)

  FCTEST_CHECK_EQUAL(connectivity%access%row(3)%col(1), 9)
  FCTEST_CHECK_EQUAL(connectivity%access%row(3)%col(2), 10)
  FCTEST_CHECK_EQUAL(connectivity%access%row(3)%col(3), 11)
  FCTEST_CHECK_EQUAL(connectivity%access%row(4)%col(1), 12)
  FCTEST_CHECK_EQUAL(connectivity%access%row(4)%col(2), 13)
  FCTEST_CHECK_EQUAL(connectivity%access%row(4)%col(3), 14)

  !============= Row access =============!

  call connectivity%row(1, row,ncols)
  FCTEST_CHECK_EQUAL(ncols,4)
  FCTEST_CHECK_EQUAL(row(1),1)
  FCTEST_CHECK_EQUAL(row(2),2)
  FCTEST_CHECK_EQUAL(row(3),3)
  FCTEST_CHECK_EQUAL(row(4),4)
  call connectivity%row(2, row,ncols)
  FCTEST_CHECK_EQUAL(ncols,4)
  FCTEST_CHECK_EQUAL(row(1),5)
  FCTEST_CHECK_EQUAL(row(2),6)
  FCTEST_CHECK_EQUAL(row(3),7)
  FCTEST_CHECK_EQUAL(row(4),8)
  call connectivity%row(3, row,ncols)
  FCTEST_CHECK_EQUAL(ncols,3)
  FCTEST_CHECK_EQUAL(row(1),9)
  FCTEST_CHECK_EQUAL(row(2),10)
  FCTEST_CHECK_EQUAL(row(3),11)
  call connectivity%row(4, row,ncols)
  FCTEST_CHECK_EQUAL(ncols,3)
  FCTEST_CHECK_EQUAL(row(1),12)
  FCTEST_CHECK_EQUAL(row(2),13)
  FCTEST_CHECK_EQUAL(row(3),14)

  call connectivity%final()
#endif
END_TEST

! -----------------------------------------------------------------------------

TEST( test_multiblockconnectivity )
  use fckit_module
  use atlas_connectivity_module
  use atlas_kinds_module
  use, intrinsic :: iso_c_binding

  implicit none
  type(atlas_MultiBlockConnectivity) :: multiblock
  type(atlas_BlockConnectivity) :: block
  integer(ATLAS_KIND_IDX), pointer :: data(:,:), padded(:,:)
  integer(ATLAS_KIND_IDX), pointer :: cols(:)

  type(atlas_Connectivity) :: base

  call fckit_log%info("test_multiblockconnectivity starting")

  multiblock = atlas_MultiBlockConnectivity()

  FCTEST_CHECK_EQUAL( multiblock%owners(), 1 )
  FCTEST_CHECK_EQUAL(multiblock%name(),"")
  FCTEST_CHECK_EQUAL(multiblock%rows(),0)
  FCTEST_CHECK_EQUAL(multiblock%blocks(),0)

  call multiblock%add(2,4, &
    & [ 1, 2, 3, 4,  &
    &   5, 6, 7, 8 ] )

  FCTEST_CHECK_EQUAL(multiblock%mincols(),4)
  FCTEST_CHECK_EQUAL(multiblock%maxcols(),4)
  FCTEST_CHECK_EQUAL(multiblock%rows(),   2)
  FCTEST_CHECK_EQUAL(multiblock%blocks(), 1)

  call multiblock%add(2,3, &
    & [ 9,  10, 11,  &
    &   12, 13, 14 ] )

  FCTEST_CHECK_EQUAL(multiblock%mincols(),3)
  FCTEST_CHECK_EQUAL(multiblock%maxcols(),4)
  FCTEST_CHECK_EQUAL(multiblock%blocks(), 2)

  block = multiblock%block(1_ATLAS_KIND_IDX)
  !FCTEST_CHECK_EQUAL( block%owners(), 2 )

  FCTEST_CHECK_EQUAL( block%rows(), 2 )
  FCTEST_CHECK_EQUAL( block%cols(), 4 )

  call block%data(data)
  FCTEST_CHECK_EQUAL(data(1,1), 1)
  FCTEST_CHECK_EQUAL(data(2,1), 2)
  FCTEST_CHECK_EQUAL(data(3,1), 3)
  FCTEST_CHECK_EQUAL(data(4,1), 4)
  FCTEST_CHECK_EQUAL(data(1,2), 5)
  FCTEST_CHECK_EQUAL(data(2,2), 6)
  FCTEST_CHECK_EQUAL(data(3,2), 7)
  FCTEST_CHECK_EQUAL(data(4,2), 8)

  block = multiblock%block(2_ATLAS_KIND_IDX)
  !FCTEST_CHECK_EQUAL( block%owners(), 2 )

  FCTEST_CHECK_EQUAL( block%rows(), 2 )
  FCTEST_CHECK_EQUAL( block%cols(), 3 )

  call block%data(data)
  FCTEST_CHECK_EQUAL(data(1,1), 9)
  FCTEST_CHECK_EQUAL(data(2,1), 10)
  FCTEST_CHECK_EQUAL(data(3,1), 11)
  FCTEST_CHECK_EQUAL(data(1,2), 12)
  FCTEST_CHECK_EQUAL(data(2,2), 13)
  FCTEST_CHECK_EQUAL(data(3,2), 14)

  call block%final()

  FCTEST_CHECK_EQUAL( multiblock%owners(), 1 )
  base = multiblock
  FCTEST_CHECK_EQUAL( base%owners(), 2 )
  FCTEST_CHECK_EQUAL( multiblock%owners(), 2 )

  call base%padded_data(padded,cols)
  FCTEST_CHECK_EQUAL(padded(1,1), 1)
  FCTEST_CHECK_EQUAL(padded(2,1), 2)
  FCTEST_CHECK_EQUAL(padded(3,1), 3)
  FCTEST_CHECK_EQUAL(padded(4,1), 4)
  FCTEST_CHECK_EQUAL(padded(1,2), 5)
  FCTEST_CHECK_EQUAL(padded(2,2), 6)
  FCTEST_CHECK_EQUAL(padded(3,2), 7)
  FCTEST_CHECK_EQUAL(padded(4,2), 8)

  FCTEST_CHECK_EQUAL(padded(1,3), 9)
  FCTEST_CHECK_EQUAL(padded(2,3), 10)
  FCTEST_CHECK_EQUAL(padded(3,3), 11)
  FCTEST_CHECK_EQUAL(padded(1,4), 12)
  FCTEST_CHECK_EQUAL(padded(2,4), 13)
  FCTEST_CHECK_EQUAL(padded(3,4), 14)

  call multiblock%final()

  FCTEST_CHECK_EQUAL( base%owners(), 1 )

  call base%final()

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

