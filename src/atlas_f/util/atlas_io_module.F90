! (C) Copyright 2013-2015 ECMWF.

#include "atlas/atlas_f.h"

module atlas_io_module

use atlas_c_interop, only : c_str
use atlas_FunctionSpace_module, only: &
    & atlas_FunctionSpace
use atlas_Field_module, only: &
    & atlas_Field
use atlas_FieldSet_module, only: &
    & atlas_FieldSet
use atlas_Mesh_module, only: &
    & atlas_Mesh

implicit none

private :: c_str
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_FieldSet
private :: atlas_Mesh

public :: atlas_read_gmsh
public :: atlas_write_gmsh
public :: atlas_write_gmsh_field
public :: atlas_write_gmsh_fieldset

public

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

! =============================================================================
CONTAINS
! =============================================================================

function atlas_read_gmsh(filename) result(mesh)
  use atlas_gmsh_c_binding
  character(len=*), intent(in) :: filename
  type(atlas_Mesh) :: mesh
  mesh = atlas_Mesh( atlas__read_gmsh(c_str(filename)) )
  call mesh%return()
end function atlas_read_gmsh

subroutine atlas_write_gmsh(mesh,filename)
  use atlas_gmsh_c_binding
  type(atlas_Mesh), intent(in) :: mesh
  character(len=*), intent(in) :: filename
  call atlas__write_gmsh_mesh(mesh%c_ptr(),c_str(filename))
end subroutine atlas_write_gmsh

subroutine atlas_write_gmsh_field(field,function_space,filename,mode)
  use atlas_gmsh_c_binding
  type(atlas_Field), intent(in) :: field
  class(atlas_functionspace), intent(in) :: function_space
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_field(field%c_ptr(),function_space%c_ptr(),c_str(filename),mode)
  else
    call atlas__write_gmsh_field(field%c_ptr(),function_space%c_ptr(),c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_field

subroutine atlas_write_gmsh_fieldset(fieldset,function_space,filename,mode)
  use atlas_gmsh_c_binding
  type(atlas_FieldSet), intent(in) :: fieldset
  class(atlas_functionspace), intent(in) :: function_space
  character(len=*), intent(in) :: filename
  integer(kind(openmode)), optional :: mode
  if( present(mode) ) then
    call atlas__write_gmsh_fieldset(fieldset%c_ptr(),function_space%c_ptr(),c_str(filename),mode)
  else
    call atlas__write_gmsh_fieldset(fieldset%c_ptr(),function_space%c_ptr(),c_str(filename),out)
  endif
end subroutine atlas_write_gmsh_fieldset

! -----------------------------------------------------------------------------

end module atlas_io_module
