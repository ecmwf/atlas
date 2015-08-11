! (C) Copyright 1996-2014 ECMWF.
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

module fctest_atlas_FunctionSpace_Fixture
use atlas_module
use iso_c_binding
implicit none

contains

end module fctest_atlas_FunctionSpace_Fixture

! -----------------------------------------------------------------------------

TESTSUITE_WITH_FIXTURE(fctest_atlas_FunctionSpace,fctest_atlas_FunctionSpace_Fixture)

! -----------------------------------------------------------------------------

TESTSUITE_INIT
  call atlas_init()
END_TESTSUITE_INIT

! -----------------------------------------------------------------------------

TESTSUITE_FINALIZE
  call atlas_finalize()
END_TESTSUITE_FINALIZE

! -----------------------------------------------------------------------------

TEST( test_nodes )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_NodesFunctionSpace) :: fs
type(atlas_Field) :: field, template
integer :: halo_size
halo_size = 1

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fs = atlas_NodesFunctionSpace(mesh,halo_size)

field = fs%create_field(atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%datatype() , "real64" )
call atlas_delete(field)

field = fs%create_field("field",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
call atlas_delete(field)

field = fs%create_field((/2/),atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_field("field",(/2,2/),atlas_integer(c_int))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call atlas_delete(field)
call atlas_delete(template)


field = fs%create_global_field(atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%datatype() , "real64" )
call atlas_delete(field)

field = fs%create_global_field("field",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 1 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
call atlas_delete(field)

field = fs%create_global_field((/2/),atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_global_field("field",(/2,2/),atlas_integer(c_int))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_global_field(template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_global_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call atlas_delete(field)
call atlas_delete(template)

call atlas_delete(fs)
call atlas_delete(mesh)
call atlas_delete(grid)

END_TEST

! -----------------------------------------------------------------------------


TEST( test_nodescolumns )
type(atlas_ReducedGrid) :: grid
type(atlas_Mesh) :: mesh
type(atlas_NodesColumnFunctionSpace) :: fs
type(atlas_Field) :: field, template
integer :: halo_size, levels
halo_size = 1
levels = 10

grid = atlas_ReducedGrid("rgg.N24")
mesh = atlas_generate_mesh(grid)
fs = atlas_NodesColumnFunctionSpace(mesh,10,halo_size)

field = fs%create_field(atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%datatype() , "real64" )
call atlas_delete(field)

field = fs%create_field("field",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
call atlas_delete(field)

field = fs%create_field((/2/),atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_field("field",(/2,2/),atlas_integer(c_int))
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_field(template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call atlas_delete(field)
call atlas_delete(template)


field = fs%create_global_field(atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "" )
FCTEST_CHECK_EQUAL( field%datatype() , "real64" )
call atlas_delete(field)

field = fs%create_global_field("field",atlas_real(c_float))
FCTEST_CHECK_EQUAL( field%rank() , 2 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
FCTEST_CHECK_EQUAL( field%datatype() , "real32" )
call atlas_delete(field)

field = fs%create_global_field((/2/),atlas_real(c_double))
FCTEST_CHECK_EQUAL( field%rank() , 3 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_global_field("field",(/2,2/),atlas_integer(c_int))
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
template = field

field = fs%create_global_field(template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "" )
call atlas_delete(field)

field = fs%create_global_field("field",template)
FCTEST_CHECK_EQUAL( field%rank() , 4 )
FCTEST_CHECK_EQUAL( field%name() , "field" )
call atlas_delete(field)
call atlas_delete(template)

call atlas_delete(fs)
call atlas_delete(mesh)
call atlas_delete(grid)

END_TEST

! -----------------------------------------------------------------------------

END_TESTSUITE

