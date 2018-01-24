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

TESTSUITE(fctest_atlas_Gmsh)

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

TEST( test_gmsh )
  use atlas_module
  implicit none
  type(atlas_StructuredGrid) :: grid
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Mesh) :: mesh
  type(atlas_functionspace_NodeColumns) :: functionspace_nodes
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_Field) :: field
  type(atlas_FieldSet) :: fieldset
  type(atlas_Output) :: gmsh
  integer(4), pointer :: fdata(:,:)
  integer :: jlev

  call atlas_log%info("test_gmsh starting")

  grid = atlas_StructuredGrid("N24")
  meshgenerator = atlas_MeshGenerator()
  mesh = meshgenerator%generate(grid)
  call meshgenerator%final()

  gmsh = atlas_output_Gmsh("output_fctest_gmsh.msh","w",coordinates="xyz",gather=.false.,levels=[0,2])
  call gmsh%write(mesh)

  functionspace_nodes = atlas_functionspace_NodeColumns(mesh,1)
  nodes = mesh%nodes()
  field = nodes%global_index()
  call gmsh%write(field,functionspace_nodes)
  field = nodes%remote_index()
  call gmsh%write(field,functionspace_nodes)
  
  fieldset = atlas_FieldSet()
  call fieldset%add( nodes%lonlat() )
  call fieldset%add( nodes%partition() )
  call gmsh%write(fieldset,functionspace_nodes)
  
  
  field = functionspace_nodes%create_field(name="leveled",kind=atlas_integer(4),levels=4);
  call field%data(fdata)
  do jlev=1,field%levels()
    fdata(jlev,:) = jlev
  enddo
  call gmsh%write(field)


  call gmsh%write(field)
  
  fieldset = atlas_FieldSet()
  field = functionspace_nodes%create_field(name="scal1",kind=atlas_integer(4),levels=4);
  call field%data(fdata)
  do jlev=1,field%levels()
    fdata(jlev,:) = jlev-1
  enddo
  call fieldset%add(field)

  field = functionspace_nodes%create_field(name="scal2",kind=atlas_integer(4),levels=4);
  call field%data(fdata)
  do jlev=1,field%levels()
    fdata(jlev,:) = -(jlev-1)
  enddo
  call fieldset%add(field)
  
  call gmsh%write(fieldset)
  
  
  END_TEST


! -----------------------------------------------------------------------------

END_TESTSUITE

