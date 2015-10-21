! (C) Copyright 1996-2015 ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! This File contains Unit Tests for testing the
! C++ / Fortran Interfaces to the State Datastructure
! @author Willem Deconinck

#include "fctest/fctest.h"

! -----------------------------------------------------------------------------

module fctest_atlas_nabla_EdgeBasedFiniteVolume_Fixture
use atlas_module
use iso_c_binding
implicit none

contains

end module fctest_atlas_nabla_EdgeBasedFiniteVolume_Fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_nabla_EdgeBasedFiniteVolume,fctest_atlas_nabla_EdgeBasedFiniteVolume_Fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_fvm )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_functionspace_EdgeBasedFiniteVolume) :: fvm

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fvm  = atlas_functionspace_EdgeBasedFiniteVolume(mesh)

call fvm%final()
call mesh%final()
call grid%final()

END_TEST

! -----------------------------------------------------------------------------

TEST( test_nabla )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_Nodes) :: nodes
type(atlas_MeshGenerator) :: meshgenerator
type(atlas_functionspace_EdgeBasedFiniteVolume) :: fvm
type(atlas_Nabla) :: nabla
type(atlas_Field) :: ghostfield
type(atlas_Field) :: varfield
type(atlas_Field) :: gradfield

real(c_double), pointer :: var(:,:)
real(c_double), pointer :: grad(:,:,:)
logical, pointer :: is_ghost(:)

integer, parameter :: nlev = 137

type(atlas_Config) :: config

config = atlas_Config()
call config%set("radius",1.0)

! Setup
grid = atlas_ReducedGrid("rgg.N24")
meshgenerator = atlas_ReducedGridMeshGenerator()
mesh = meshgenerator%generate(grid) ! second optional argument for atlas_GridDistrubution
fvm  = atlas_functionspace_EdgeBasedFiniteVolume(mesh,config)
nabla = atlas_Nabla(fvm)

! Create a variable field and a gradient field
varfield = fvm%create_field("var",atlas_real(c_double),nlev)
gradfield  = fvm%create_field("grad",atlas_real(c_double),nlev,[2])

! Access to data
call varfield%data(var)
call gradfield%data(grad)
var(:,:) = 0.

call fvm%halo_exchange(varfield)

! Compute the gradient
call nabla%gradient(varfield,gradfield)

! get is_ghost
nodes = mesh%nodes()
ghostfield = nodes%ghost()
call ghostfield%data(is_ghost)

! Cleanup
call config%final()
call ghostfield%final()
call varfield%final()
call gradfield%final()
call nabla%final()
call fvm%final()
call nodes%final()
call mesh%final()
call grid%final()
call meshgenerator%final()

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

