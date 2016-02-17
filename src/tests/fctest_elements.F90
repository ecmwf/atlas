! (C) Copyright 1996-2016 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the Elements
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

TESTSUITE(fctest_atlas_elements)

! -----------------------------------------------------------------------------

TEST( test_elementtype )
  use atlas_mesh_hybridelements_module
  use atlas_mesh_ElementType_module
  use atlas_mesh_Elements_module

  implicit none
  type(atlas_mesh_ElementType) :: triag, quad, line

  write(*,*) "test_elementtype"

  triag = atlas_mesh_Triangle()
  quad  = atlas_mesh_Quadrilateral()
  line  = atlas_mesh_Line()

END_TEST

! -----------------------------------------------------------------------------

TEST( test_hybridelements )
  use atlas_mesh_hybridelements_module
  use atlas_mesh_ElementType_module
  use atlas_mesh_Elements_module
  use atlas_connectivity_module
  use atlas_field_module
  use, intrinsic :: iso_c_binding

  implicit none
  type(atlas_mesh_Cells) :: cells
  type(atlas_MultiBlockConnectivity) :: node_connectivity
  type(atlas_Field) :: field, field2
  type(atlas_mesh_Elements) :: elements
  type(atlas_mesh_ElementType) :: element_type
  integer(c_int), pointer :: data(:,:)
  integer(c_size_t) :: jfield

  write(*,*) "test_hybridelements starting"

  cells = atlas_mesh_Cells()
  FCTEST_CHECK_EQUAL( cells%owners(), 1 )

  FCTEST_CHECK_EQUAL( cells%size(), 0_c_size_t )
  node_connectivity = cells%node_connectivity()
  FCTEST_CHECK_EQUAL( node_connectivity%owners(), 2 )
  FCTEST_CHECK_EQUAL( node_connectivity%rows(), 0_c_size_t )
  FCTEST_CHECK( cells%has_field("glb_idx") )
  FCTEST_CHECK( cells%has_field("partition") )
  FCTEST_CHECK( cells%has_field("remote_idx") )
  FCTEST_CHECK( cells%has_field("halo") )

  FCTEST_CHECK_EQUAL( cells%nb_fields(), 4_c_size_t )
  field = cells%field("partition")
  FCTEST_CHECK_EQUAL( field%owners(), 2 )
  FCTEST_CHECK_EQUAL( field%name(), "partition" )

  field2 = cells%partition()
  FCTEST_CHECK_EQUAL( field2%name(), "partition" )
  FCTEST_CHECK_EQUAL( field%owners(), 3 )

  field2 = cells%global_index()
  FCTEST_CHECK_EQUAL( field2%name(), "glb_idx" )
  FCTEST_CHECK_EQUAL( field%owners(), 2 )

  field2 = cells%remote_index()
  FCTEST_CHECK_EQUAL( field2%name(), "remote_idx" )
  FCTEST_CHECK_EQUAL( field%owners(), 2 )

  field2 = cells%halo()
  FCTEST_CHECK_EQUAL( field2%name(), "halo" )
  FCTEST_CHECK_EQUAL( field%owners(), 2 )

  call field2%final()
  do jfield=1,cells%nb_fields()
    field = cells%field(jfield)
    write(0,*) field%name()
    FCTEST_CHECK_EQUAL( field%owners(), 2 )
  enddo



  call cells%add( atlas_mesh_Triangle(), 5_c_size_t , &
        &  [  1,  2  ,3, &
        &     4,  5,  6, &
        &     7,  8,  9, &
        &    10, 11, 12, &
        &    13, 14, 15  ] )


  FCTEST_CHECK_EQUAL( cells%size(), 5_c_size_t )
  FCTEST_CHECK_EQUAL( cells%nb_types(), 1_c_size_t )
  FCTEST_CHECK_EQUAL( node_connectivity%rows(),    5_c_size_t )
  FCTEST_CHECK_EQUAL( node_connectivity%maxcols(), 3_c_size_t )

  call node_connectivity%data(data)

  FCTEST_CHECK_EQUAL( data(1,1) ,  1 )
  FCTEST_CHECK_EQUAL( data(2,1) ,  2 )
  FCTEST_CHECK_EQUAL( data(3,1) ,  3 )
  FCTEST_CHECK_EQUAL( data(1,2) ,  4 )
  FCTEST_CHECK_EQUAL( data(2,2) ,  5 )
  FCTEST_CHECK_EQUAL( data(3,2) ,  6 )
  FCTEST_CHECK_EQUAL( data(1,3) ,  7 )
  FCTEST_CHECK_EQUAL( data(2,3) ,  8 )
  FCTEST_CHECK_EQUAL( data(3,3) ,  9 )
  FCTEST_CHECK_EQUAL( data(1,4) , 10 )
  FCTEST_CHECK_EQUAL( data(2,4) , 11 )
  FCTEST_CHECK_EQUAL( data(3,4) , 12 )
  FCTEST_CHECK_EQUAL( data(1,5) , 13 )
  FCTEST_CHECK_EQUAL( data(2,5) , 14 )
  FCTEST_CHECK_EQUAL( data(3,5) , 15 )

  call cells%add( atlas_mesh_Quadrilateral(), 2_c_size_t , &
        &  [  16, 17, 18, 19, &
        &     20, 21, 22, 23  ] )

  FCTEST_CHECK_EQUAL( cells%nb_types(), 2_c_size_t )
  FCTEST_CHECK_EQUAL( cells%size(), 7_c_size_t )
  FCTEST_CHECK_EQUAL( node_connectivity%rows(),    7_c_size_t )
  FCTEST_CHECK_EQUAL( node_connectivity%mincols(), 3_c_size_t )
  FCTEST_CHECK_EQUAL( node_connectivity%maxcols(), 4_c_size_t )

  elements = cells%elements(1)
  FCTEST_CHECK_EQUAL( elements%owners(), 2 )
  FCTEST_CHECK_EQUAL( elements%size(), 5_c_size_t )
  element_type = elements%element_type()

  ! Should print ERROR
  call node_connectivity%data(data)

  call node_connectivity%padded_data(data)
  FCTEST_CHECK_EQUAL( data(1,6) , 16 )
  FCTEST_CHECK_EQUAL( data(2,6) , 17 )
  FCTEST_CHECK_EQUAL( data(3,6) , 18 )
  FCTEST_CHECK_EQUAL( data(4,6) , 19 )
  FCTEST_CHECK_EQUAL( data(1,7) , 20 )
  FCTEST_CHECK_EQUAL( data(2,7) , 21 )
  FCTEST_CHECK_EQUAL( data(3,7) , 22 )
  FCTEST_CHECK_EQUAL( data(4,7) , 23 )


  field = atlas_Field("p",atlas_real(8),[int(cells%size())])
  call cells%add(field)

  call cells%final()

  FCTEST_CHECK_EQUAL( node_connectivity%owners(), 1 )
  FCTEST_CHECK_EQUAL( field%owners(), 1 )


  call node_connectivity%final()
  call field%final()

END_TEST

! -----------------------------------------------------------------------------

TEST( test_elements )
  use atlas_mesh_hybridelements_module
  use atlas_mesh_ElementType_module
  use atlas_mesh_Elements_module
  use atlas_connectivity_module
  use atlas_field_module
  use, intrinsic :: iso_c_binding

  implicit none
  type(atlas_mesh_Cells) :: cells
  type(atlas_BlockConnectivity) :: node_connectivity
  type(atlas_mesh_Elements) :: elements
  type(atlas_mesh_ElementType) :: element_type
  integer(c_int), pointer :: data(:,:)

  write(*,*) "test_elements starting"

  cells = atlas_mesh_Cells()
  call cells%add( atlas_mesh_Triangle(), 5_c_size_t , &
        &  [  1,  2  ,3, &
        &     4,  5,  6, &
        &     7,  8,  9, &
        &    10, 11, 12, &
        &    13, 14, 15  ] )

  call cells%add( atlas_mesh_Quadrilateral(), 2_c_size_t , &
        &  [  16, 17, 18, 19, &
        &     20, 21, 22, 23  ] )

  elements = cells%elements(1)

  FCTEST_CHECK_EQUAL( elements%begin(), 1_c_size_t )
  FCTEST_CHECK_EQUAL( elements%end(), 5_c_size_t )

  element_type = elements%element_type()
  FCTEST_CHECK_EQUAL( element_type%owners(), 2 )

  node_connectivity = elements%node_connectivity()
  FCTEST_CHECK_EQUAL( node_connectivity%owners(), 2 )

  FCTEST_CHECK_EQUAL( element_type%nb_nodes(), 3_c_size_t )
  FCTEST_CHECK_EQUAL( element_type%nb_edges(), 3_c_size_t )
  FCTEST_CHECK_EQUAL( element_type%name(), "Triangle" )
  FCTEST_CHECK( element_type%parametric() )


  call node_connectivity%data(data)

  FCTEST_CHECK_EQUAL( data(1,1) ,  1 )
  FCTEST_CHECK_EQUAL( data(2,1) ,  2 )
  FCTEST_CHECK_EQUAL( data(3,1) ,  3 )
  FCTEST_CHECK_EQUAL( data(1,2) ,  4 )
  FCTEST_CHECK_EQUAL( data(2,2) ,  5 )
  FCTEST_CHECK_EQUAL( data(3,2) ,  6 )
  FCTEST_CHECK_EQUAL( data(1,3) ,  7 )
  FCTEST_CHECK_EQUAL( data(2,3) ,  8 )
  FCTEST_CHECK_EQUAL( data(3,3) ,  9 )
  FCTEST_CHECK_EQUAL( data(1,4) , 10 )
  FCTEST_CHECK_EQUAL( data(2,4) , 11 )
  FCTEST_CHECK_EQUAL( data(3,4) , 12 )
  FCTEST_CHECK_EQUAL( data(1,5) , 13 )
  FCTEST_CHECK_EQUAL( data(2,5) , 14 )
  FCTEST_CHECK_EQUAL( data(3,5) , 15 )

  elements = cells%elements(2)

  FCTEST_CHECK_EQUAL( elements%begin(), 6_c_size_t )
  FCTEST_CHECK_EQUAL( elements%end(), 7_c_size_t )

  element_type = elements%element_type()
  FCTEST_CHECK_EQUAL( element_type%owners(), 2 )

  FCTEST_CHECK_EQUAL( element_type%nb_nodes(), 4_c_size_t )
  FCTEST_CHECK_EQUAL( element_type%nb_edges(), 4_c_size_t )
  FCTEST_CHECK_EQUAL( element_type%name(), "Quadrilateral" )
  FCTEST_CHECK( element_type%parametric() )

  node_connectivity = elements%node_connectivity()
  FCTEST_CHECK_EQUAL( node_connectivity%owners(), 2 )

  call node_connectivity%data(data)

  FCTEST_CHECK_EQUAL( data(1,1) , 16 )
  FCTEST_CHECK_EQUAL( data(2,1) , 17 )
  FCTEST_CHECK_EQUAL( data(3,1) , 18 )
  FCTEST_CHECK_EQUAL( data(4,1) , 19 )
  FCTEST_CHECK_EQUAL( data(1,2) , 20 )
  FCTEST_CHECK_EQUAL( data(2,2) , 21 )
  FCTEST_CHECK_EQUAL( data(3,2) , 22 )
  FCTEST_CHECK_EQUAL( data(4,2) , 23 )

!--------------------------


  ! Add elements to triangles
  elements = cells%elements(1)
  call elements%add(2)
  FCTEST_CHECK_EQUAL( elements%begin(), 1_c_size_t )
  FCTEST_CHECK_EQUAL( elements%end(), 7_c_size_t )

  node_connectivity = elements%node_connectivity()
  call node_connectivity%data(data)

  FCTEST_CHECK_EQUAL( data(1,6) , node_connectivity%missing_value() )
  FCTEST_CHECK_EQUAL( data(2,6) , node_connectivity%missing_value() )
  FCTEST_CHECK_EQUAL( data(3,6) , node_connectivity%missing_value() )
  FCTEST_CHECK_EQUAL( data(1,7) , node_connectivity%missing_value() )
  FCTEST_CHECK_EQUAL( data(2,7) , node_connectivity%missing_value() )
  FCTEST_CHECK_EQUAL( data(3,7) , node_connectivity%missing_value() )


  elements = cells%elements(2)
  FCTEST_CHECK_EQUAL( elements%begin(), 8_c_size_t )
  FCTEST_CHECK_EQUAL( elements%end(), 9_c_size_t )


  call elements%final()
  call cells%final()
  call node_connectivity%final()
  call element_type%final()

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

