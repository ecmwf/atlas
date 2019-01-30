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

TESTSUITE(fctest_atlas_Meshgen)

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

!DISABLE( from_constructor )

!  use atlas_module
!  implicit none

!  type(atlas_StructuredGrid) :: grid
!  type(atlas_Mesh) :: mesh
!  type(atlas_MeshGenerator) :: meshgenerator
!   type(atlas_Output) :: gmsh
!  integer, parameter :: nlat = 7
!  real(c_double) :: lats(nlat)
!  real(c_double) :: lon_min(nlat)
!  integer(c_int) :: pl(nlat)
!  lats    = [ 90., 60., 30.,   0., -30., -60., -90. ]
!  pl      = [  2,  4,  8,  12,   8,   4,   2 ]
!  lon_min = [  0., 15.,  0.,  30.,   0.,  15.,   0. ]
!  grid = atlas_grid_CustomStructured(lats,pl,lon_min)
!  meshgenerator = atlas_MeshGenerator()
!  mesh = meshgenerator%generate(grid)
! gmsh = atlas_output_Gmsh("test_custom_structured1.msh")
! call gmsh%write(mesh)
!  call meshgenerator%final()
!  call grid%final()
!  call mesh%final()

!END_DISABLE


!DISABLE( from_config )

!  use atlas_module
!  implicit none

!  type(atlas_Config) :: conf
!  type(atlas_StructuredGrid) :: grid
!  type(atlas_Mesh) :: mesh
!  type(atlas_MeshGenerator) :: meshgenerator
!   type(atlas_Output) :: gmsh
!  conf = atlas_Config( )
!  call conf%set("Implementationype","custom_structured")
!  call conf%set("nlat",7)
!  call conf%set("latitudes",[ 90., 60., 30.,   0., -30., -60., -90. ])
!  call conf%set("pl",[  2,  4,  8,  12,   8,   4,   2 ])
!  call conf%set("lon_min",[  0., 15.,  0.,  30.,   0.,  15.,   0. ])

!  grid = atlas_StructuredGrid(conf)
!  meshgenerator = atlas_MeshGenerator()
!  mesh = meshgenerator%generate(grid)
! gmsh = atlas_output_Gmsh("test_custom_structured2.msh")
! call gmsh%write(mesh)
!  call meshgenerator%final()
!  call grid%final()
!  call mesh%final()
!  call conf%final()
!END_DISABLE


!DISABLE( from_json_file )

!  use atlas_module
!  implicit none

!  type(atlas_Config) :: conf
!  type(atlas_StructuredGrid) :: grid
!  type(atlas_Mesh) :: mesh
!  type(atlas_MeshGenerator) :: meshgenerator
!   type(atlas_Output) :: gmsh

! ! Write a json file
! OPEN (UNIT=9 , FILE="custom.json", STATUS='REPLACE')
! write(9,'(A)') &
!   &     '{' &
!   & //  '  "Implementationype" : "custom_structured",'               &
!   & //  '  "nlat"      : 7,'                                 &
!   & //  '  "latitudes" : [ 90, 60, 30,   0, -30, -60, -90 ],' &
!   & //  '  "pl"        : [  2,  4,  8,  12,   8,   4,   2 ],' &
!   & //  '  "lon_min"   : [  0, 15,  0,  30,   0,  15,   0 ]'  &
!   & //  '}'
! CLOSE(9)

!  conf = atlas_Config( atlas_PathName("custom.json") )
!  grid = atlas_StructuredGrid(conf)
!  meshgenerator = atlas_MeshGenerator()
!  mesh = meshgenerator%generate(grid)
! gmsh = atlas_output_Gmsh("test_custom_structured3.msh")
! call gmsh%write(mesh)
!  call meshgenerator%final()
!  call grid%final()
!  call mesh%final()
!  call conf%final()
!END_DISABLE


TEST( test_meshgen )
  use, intrinsic :: iso_c_binding
  use atlas_module
  implicit none
  type(atlas_StructuredGrid) :: grid
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Mesh) :: mesh
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_mesh_Edges) :: edges
  type(atlas_Field) :: field
  type(atlas_functionspace_NodeColumns) :: functionspace_nodes
  type(atlas_functionspace_EdgeColumns) :: functionspace_edges
  type(atlas_HaloExchange) :: halo_exchange
  type(atlas_Output) :: gmsh
  integer(ATLAS_KIND_IDX), pointer :: ridx(:)
  real(c_double), pointer :: arr1d(:), arr2d(:,:)
  integer :: i, nnodes, nghost

  write(*,*) "test_meshgen starting"

  grid = atlas_StructuredGrid("N24")
!  grid = atlas_grid_ShiftedLat(40,20)

  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid)
  nodes = mesh%nodes()
  call meshgenerator%final()

  functionspace_edges = atlas_functionspace_EdgeColumns(mesh,halo=1)
  functionspace_nodes = atlas_functionspace_NodeColumns(mesh,halo=1)
  call atlas_build_median_dual_mesh(mesh)

  nnodes = nodes%size()

  FCTEST_CHECK_EQUAL( nnodes, 3672 )

  field = nodes%field("remote_idx")
  call field%data(ridx)
  nghost = 0
  do i=1,nnodes
    if( ridx(i) /= i ) nghost = nghost + 1
  enddo

  write(0,*) "nghost =",nghost

  FCTEST_CHECK_EQUAL( nghost, 144 )


  call atlas_log%info( nodes%str() )

  field = nodes%field("dual_volumes")
  call field%data(arr1d)
  call field%final()

  functionspace_nodes = atlas_functionspace_NodeColumns(mesh,1)
  halo_exchange = functionspace_nodes%get_halo_exchange()
  call halo_exchange%execute(arr1d)


  edges = mesh%edges()
  field = edges%field("dual_normals")
  call field%data(arr2d)
  call field%final()
  gmsh = atlas_output_Gmsh("testf2.msh",ghost=.true.)
  call gmsh%write(mesh)

  call atlas_write_load_balance_report(mesh,"N24_loadbalance.dat")
END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

