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

TESTSUITE(fctest_atlas_TriangularMeshBuilder)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  use atlas_module
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  use atlas_module
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_mesbuilder )
  use, intrinsic :: iso_c_binding
  use atlas_module
  implicit none
  type(atlas_Mesh) :: mesh

  integer :: nb_nodes
  integer :: nb_triags
  integer(ATLAS_KIND_GIDX), allocatable :: node_global_index(:)
  integer(ATLAS_KIND_GIDX), allocatable :: triag_global_index(:)
  integer(ATLAS_KIND_GIDX), allocatable :: triag_nodes_global_index(:,:)
  real(c_double), allocatable :: xy_in(:,:), lonlat_in(:,:)

  integer(ATLAS_KIND_GIDX) :: global_index_base = 0

  !   small mesh
  !
  !      1 ---- 5 ------- 6
  !      | 3 / 4 \ 1  / 2 |
  !      2 ------ 3 ----- 4
  !
  nb_nodes = 6
  nb_triags = 4

  allocate(node_global_index(nb_nodes))
  allocate(xy_in(2,nb_nodes))
  allocate(lonlat_in(2,nb_nodes))
  allocate(triag_global_index(nb_triags))
  allocate(triag_nodes_global_index(3,nb_triags))

  node_global_index = [6, 5, 4, 3, 2, 1]
  lonlat_in = reshape([&
    15.0, 5.0,& ! point 6
     5.0, 5.0,& ! point 5
    15.0, 0.0,& ! point 4
    10.0, 0.0,& ! point 3
     0.0, 0.0,& ! point 2
     0.0, 5.0 & ! point 1
  ], shape(lonlat_in))

  xy_in(:,:) = lonlat_in / 10._c_double

  triag_global_index = [1, 2, 3, 4]
  triag_nodes_global_index = reshape([3,6,5 , 3,4,6,  2,5,1,  2,3,5], shape(triag_nodes_global_index))

  ! Change base
  node_global_index = node_global_index - 1 + global_index_base
  triag_global_index = triag_global_index - 1 + global_index_base
  triag_nodes_global_index = triag_nodes_global_index - 1 + global_index_base

  ! Manually build mesh
  BLOCK
    type(atlas_TriangularMeshBuilder) :: meshbuilder
    meshbuilder = atlas_TriangularMeshBuilder()
    mesh = meshbuilder%build(nb_nodes, node_global_index, xy_in(1,:), xy_in(2,:), lonlat_in(1,:), lonlat_in(2,:), &
                             nb_triags, triag_global_index, triag_nodes_global_index, &
                             global_index_base)
    call meshbuilder%final()
  END BLOCK

  ! Mesh should created now. Verify some fields
  BLOCK
    type(atlas_mesh_Nodes) :: nodes
    real(c_double), pointer :: xy(:,:)
    real(c_double), pointer :: lonlat(:,:)
    integer(ATLAS_KIND_GIDX), pointer :: global_index(:)

    type(atlas_Field) :: field_xy
    type(atlas_Field) :: field_lonlat
    type(atlas_Field) :: field_global_index
 
    integer :: jnode
    integer :: nb_nodes
    nodes = mesh%nodes()
    nb_nodes = nodes%size()
    FCTEST_CHECK_EQUAL( nb_nodes, 6 )
    field_xy = nodes%xy()
    field_lonlat = nodes%lonlat()
    field_global_index = nodes%global_index()
    call field_xy%data(xy)
    call field_lonlat%data(lonlat)
    call field_global_index%data(global_index)
    do jnode=1,nb_nodes
      FCTEST_CHECK_EQUAL( xy(:,jnode), xy_in(:,jnode) )
      FCTEST_CHECK_EQUAL( lonlat(:,jnode), lonlat_in(:,jnode) )
      FCTEST_CHECK_EQUAL( global_index(jnode), node_global_index(jnode) - global_index_base + 1 )
    enddo
    call field_xy%final()
    call field_lonlat%final()
    call field_global_index%final()
    call nodes%final()
  END BLOCK

  ! Output mesh
  BLOCK
    type(atlas_Output) :: gmsh
    gmsh = atlas_output_Gmsh("out.msh",coordinates="lonlat")
    call gmsh%write(mesh)
    call gmsh%final()
  END BLOCK


  call mesh%final()

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

