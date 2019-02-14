! (C) Copyright 2013 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Mesh Datastructure
! @author Willem Deconinck

#include "fckit/fctest.h"

! -----------------------------------------------------------------------------

module fcta_Mesh_fixture
use atlas_module
use, intrinsic :: iso_c_binding
implicit none
end module fcta_Mesh_fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_Mesh,fcta_Mesh_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

!!! WARNING ! THIS IS DEPRECATED AND SHOULD NOT BE USED AS EXAMPLE !!!!

TEST( test_mesh_nodes )
implicit none

  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  integer(c_int) :: nb_nodes

  write(*,*) "test_function_space starting"
  mesh = atlas_Mesh()
  nodes = mesh%nodes()
  nb_nodes = nodes%size()
  FCTEST_CHECK_EQUAL( nb_nodes, 0 )
  FCTEST_CHECK_EQUAL( nodes%size() , 0  )
  FCTEST_CHECK( nodes%has_field("partition") )
  FCTEST_CHECK( nodes%has_field("remote_idx") )
  call nodes%resize(10)
  nb_nodes = nodes%size()
  FCTEST_CHECK_EQUAL( nb_nodes, 10 )
  FCTEST_CHECK_EQUAL( nodes%size() , 10  )
  call atlas_log%info( nodes%str() )

  call mesh%final()
  call nodes%final()
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

