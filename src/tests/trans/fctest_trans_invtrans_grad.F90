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

module fctest_atlas_trans_igrad_fixture
use atlas_module
use iso_c_binding
implicit none
character(len=1024) :: msg
real(c_double), parameter :: pi = 2.0_c_double*asin(1.0_c_double)
real(c_double), parameter :: deg2rad = pi/180._c_double
contains

! @brief Compute magnitude of flow with rotation-angle beta
! (beta=0 --> zonal, beta=pi/2 --> meridional)
subroutine rotated_flow_magnitude(fs,field,beta)
  use, intrinsic :: iso_c_binding
  type(atlas_functionspace_NodeColumns) :: fs
  type(atlas_Field) :: field
  real(c_double) :: beta

  real(c_double) :: radius
  real(c_double) :: USCAL
  real(c_double) :: pvel
  real(c_double) :: pi
  real(c_double) :: deg2rad

  real(c_double), pointer :: lonlat_deg(:,:)
  real(c_double), pointer :: var(:)
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_Field) :: lonlat
  integer :: jnode
  real(c_double) :: x,y, Ux, Uy

  radius = 6371229._c_double
  USCAL  = 20._c_double
  pvel = USCAL/radius
  pi = 2.0_c_double*asin(1.0_c_double)
  deg2rad = pi/180._c_double

  nodes = fs%nodes()
  lonlat = nodes%lonlat()
  call lonlat%data(lonlat_deg)
  call field%data(var)

  do jnode = 1, nodes%size()
    x = lonlat_deg(1,jnode) * deg2rad
    y = lonlat_deg(2,jnode) * deg2rad
    Ux =  pvel*(cos(beta)+tan(y)*cos(x)*sin(beta))*radius*cos(y)
    Uy = -pvel*sin(x)*sin(beta)*radius
    var(jnode) = sqrt(Ux*Ux+Uy*Uy)
  enddo

  call lonlat%final()
  call nodes%final()

end subroutine

end module

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_trans_igrad,fctest_atlas_trans_igrad_fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_library%initialise()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  write(0,*) "FINALIZE"
  call atlas_library%finalise()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_trans_invtrans_grad )
  type(atlas_GaussianGrid) :: grid
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Mesh) :: mesh
  type(atlas_Trans) :: trans
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_functionspace_NodeColumns) :: nodes_fs
  type(atlas_functionspace_Spectral) :: spectral_fs
  type(atlas_Field)         :: scalar
  type(atlas_Field)         :: scalar_sp
  type(atlas_Field)         :: grad
  type(atlas_Output)        :: gmsh

  real(c_double) :: beta

  grid = atlas_GaussianGrid("O32")
  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid)
  call meshgenerator%final()

  trans = atlas_Trans(grid, int(grid%N()-1) )

  nodes_fs    = atlas_functionspace_NodeColumns(mesh,0)
  spectral_fs = atlas_functionspace_Spectral(trans%truncation())

  scalar = nodes_fs%create_field(name="scalar",kind=atlas_real(c_double))
  grad   = nodes_fs%create_field(name="grad",kind=atlas_real(c_double),variables=2)
  scalar_sp = spectral_fs%create_field(name="spectral_sp",kind=atlas_real(c_double))

  beta = pi*0.5;
  call rotated_flow_magnitude(nodes_fs,scalar,beta)

  call trans%dirtrans(scalar,scalar_sp)

  call trans%invtrans_grad(scalar_sp,grad);

  call nodes_fs%halo_exchange(grad);


! output
  gmsh = atlas_output_Gmsh("O32-nodes-fort.msh")
  call gmsh%write(mesh)
  call gmsh%write(scalar)
  call gmsh%write(grad)

! cleanup
  call gmsh%final()

  call nodes_fs%final()
  call spectral_fs%final()

  call scalar%final()
  call scalar_sp%final()
  call grad%final()
  call nodes%final()
  call mesh%final()
  call trans%final()
  call grid%final()
END_TEST

END_TESTSUITE

