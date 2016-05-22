! (C) Copyright 2013-2015 ECMWF.

#include "atlas/atlas_f.h"

module atlas_io_module

use atlas_c_interop, only : c_str
use atlas_refcounted_module, only : atlas_refcounted
use atlas_Config_module, only : atlas_Config
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
private :: atlas_Config
private :: atlas_refcounted

public :: atlas_read_gmsh
public :: atlas_write_gmsh
public :: atlas_write_gmsh_field
public :: atlas_write_gmsh_fieldset

public :: atlas_Output
public :: atlas_output_Gmsh

public

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM



#if 1
!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_Output

! Purpose :
! -------

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_Output__delete
  procedure, public :: copy => atlas_Output__copy
  procedure, private :: write_mesh
  procedure, private :: write_field_fs
  procedure, private :: write_field
  procedure, private :: write_fieldset_fs
  procedure, private :: write_fieldset
  generic, public :: write => &
    & write_mesh, &
    & write_field_fs, &
    & write_fieldset_fs, &
    & write_field, &
    & write_fieldset

END TYPE atlas_Output

interface atlas_Output
  module procedure atlas_Output__cptr
end interface

interface atlas_output_Gmsh
  module procedure atlas_output_Gmsh__pathname_mode
end interface

#endif

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


function atlas_Output__cptr(cptr) result(Output)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_Output) :: Output
  type(c_ptr), intent(in) :: cptr
  call Output%reset_c_ptr( cptr )
end function

subroutine atlas_Output__delete(this)
  use atlas_Output_c_binding
  class(atlas_Output), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__Output__delete(this%c_ptr())
  endif
  call this%reset_c_ptr()
end subroutine atlas_Output__delete


subroutine atlas_Output__copy(this,obj_in)
  class(atlas_Output), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine

function atlas_output_Gmsh__pathname_mode(file,mode,coordinates,levels,gather) result(output)
  use atlas_output_gmsh_c_binding
  type(atlas_Output) :: output
  character(len=*), intent(in) :: file
  character(len=1), intent(in), optional :: mode
  character(len=*), intent(in), optional :: coordinates
  integer, intent(in), optional :: levels(:)
  logical, intent(in), optional :: gather
  character(len=1) :: opt_mode
  type(atlas_Config) :: opt_config
  opt_config = atlas_Config()
  opt_mode = "w"
  if( present(mode) ) opt_mode = mode
  if( present(coordinates) ) call opt_config%set("coordinates",coordinates)
  if( present(levels) )      call opt_config%set("levels",levels)
  if( present(gather) )      call opt_config%set("gather",gather)
  output = atlas_Output__cptr(atlas__output__Gmsh__create_pathname_mode_config(c_str(file),c_str(opt_mode),opt_config%c_ptr()))
  call output%return()
  call opt_config%final()
end function

subroutine write_mesh(this,mesh,config)
  use atlas_output_c_binding
  class(atlas_Output), intent(in) :: this
  type(atlas_Mesh), intent(in) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call atlas__output__write_mesh(this%c_ptr(),mesh%c_ptr(),config%c_ptr())
  else
    opt_config = atlas_Config()
    call atlas__output__write_mesh(this%c_ptr(),mesh%c_ptr(),opt_config%c_ptr())
    call opt_config%final()
  endif
end subroutine

subroutine write_field_fs(this,field,functionspace,config)
  use atlas_output_c_binding
  class(atlas_Output), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  class(atlas_FunctionSpace), intent(in) :: functionspace
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call atlas__output__write_field_fs(this%c_ptr(),field%c_ptr(),functionspace%c_ptr(),config%c_ptr())
  else
    opt_config = atlas_Config()
    call atlas__output__write_field_fs(this%c_ptr(),field%c_ptr(),functionspace%c_ptr(),opt_config%c_ptr())
    call opt_config%final()
  endif
end subroutine

subroutine write_fieldset_fs(this,fieldset,functionspace,config)
  use atlas_output_c_binding
  class(atlas_Output), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  class(atlas_FunctionSpace), intent(in) :: functionspace
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call atlas__output__write_fieldset_fs(this%c_ptr(),fieldset%c_ptr(),functionspace%c_ptr(),config%c_ptr())
  else
    opt_config = atlas_Config()
    call atlas__output__write_fieldset_fs(this%c_ptr(),fieldset%c_ptr(),functionspace%c_ptr(),opt_config%c_ptr())
    call opt_config%final()
  endif
end subroutine


subroutine write_field(this,field,config)
  use atlas_output_c_binding
  class(atlas_Output), intent(in) :: this
  type(atlas_Field), intent(in) :: field
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call atlas__output__write_field(this%c_ptr(),field%c_ptr(),config%c_ptr())
  else
    opt_config = atlas_Config()
    call atlas__output__write_field(this%c_ptr(),field%c_ptr(),opt_config%c_ptr())
    call opt_config%final()
  endif
end subroutine


subroutine write_fieldset(this,fieldset,config)
  use atlas_output_c_binding
  class(atlas_Output), intent(in) :: this
  type(atlas_FieldSet), intent(in) :: fieldset
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call atlas__output__write_fieldset(this%c_ptr(),fieldset%c_ptr(),config%c_ptr())
  else
    opt_config = atlas_Config()
    call atlas__output__write_fieldset(this%c_ptr(),fieldset%c_ptr(),opt_config%c_ptr())
    call opt_config%final()
  endif
end subroutine

end module atlas_io_module
