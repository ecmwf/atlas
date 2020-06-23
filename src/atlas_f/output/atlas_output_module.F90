! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_output_module

use fckit_owned_object_module, only : fckit_owned_object
use atlas_Config_module, only : atlas_Config
use atlas_FunctionSpace_module, only: atlas_FunctionSpace
use atlas_FieldSet_module, only: atlas_FieldSet
use atlas_Field_module, only: atlas_Field
use atlas_Mesh_module, only: atlas_Mesh
  use, intrinsic :: iso_c_binding, only : c_ptr

implicit none

public :: atlas_Output
public :: atlas_output_Gmsh

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Output

! Purpose :
! -------

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
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

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Output__final_auto
#endif

END TYPE atlas_Output

interface atlas_Output
  module procedure atlas_Output__cptr
end interface

interface atlas_output_Gmsh
  module procedure atlas_output_Gmsh__pathname_mode
end interface

!------------------------------------------------------------------------------

private :: fckit_owned_object
private :: atlas_Config
private :: atlas_FunctionSpace
private :: atlas_FieldSet
private :: atlas_Field
private :: atlas_Mesh
private :: c_ptr

! =============================================================================
CONTAINS
! =============================================================================

function atlas_Output__cptr(cptr) result(this)
  type(atlas_Output) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_output_Gmsh__pathname_mode(file,mode,coordinates,levels,gather,ghost) result(this)
  use fckit_c_interop_module, only : c_str
  use atlas_output_gmsh_c_binding
  type(atlas_Output) :: this
  character(len=*), intent(in) :: file
  character(len=1), intent(in), optional :: mode
  character(len=*), intent(in), optional :: coordinates
  integer, intent(in), optional :: levels(:)
  logical, intent(in), optional :: gather
  logical, intent(in), optional :: ghost
  character(len=1) :: opt_mode
  type(atlas_Config) :: opt_config
  opt_config = atlas_Config()
  opt_mode = "w"
  if( present(mode) ) opt_mode = mode
  if( present(coordinates) ) call opt_config%set("coordinates",coordinates)
  if( present(levels) )      call opt_config%set("levels",levels)
  if( present(gather) )      call opt_config%set("gather",gather)
  if( present(ghost) )       call opt_config%set("ghost",ghost)
  call this%reset_c_ptr( atlas__output__Gmsh__create_pathname_mode_config(c_str(file),c_str(opt_mode),&
     opt_config%CPTR_PGIBUG_B) )
  call this%return()
  call opt_config%final()
end function

subroutine write_mesh(this,mesh,config)
  use atlas_output_c_binding
  class(atlas_Output), intent(in) :: this
  type(atlas_Mesh), intent(in) :: mesh
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    call atlas__output__write_mesh(this%CPTR_PGIBUG_A,mesh%CPTR_PGIBUG_A,config%CPTR_PGIBUG_B)
  else
    opt_config = atlas_Config()
    call atlas__output__write_mesh(this%CPTR_PGIBUG_A,mesh%CPTR_PGIBUG_A,opt_config%CPTR_PGIBUG_B)
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
    call atlas__output__write_field_fs(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A,functionspace%CPTR_PGIBUG_A, &
       config%CPTR_PGIBUG_B)
  else
    opt_config = atlas_Config()
    call atlas__output__write_field_fs(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A,functionspace%CPTR_PGIBUG_A, &
      opt_config%CPTR_PGIBUG_B)
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
    call atlas__output__write_fieldset_fs(this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A,functionspace%CPTR_PGIBUG_A, &
      config%CPTR_PGIBUG_B)
  else
    opt_config = atlas_Config()
    call atlas__output__write_fieldset_fs(this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A,functionspace%CPTR_PGIBUG_A, &
      opt_config%CPTR_PGIBUG_B)
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
    call atlas__output__write_field(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A, &
      config%CPTR_PGIBUG_B)
  else
    opt_config = atlas_Config()
    call atlas__output__write_field(this%CPTR_PGIBUG_A,field%CPTR_PGIBUG_A, &
      opt_config%CPTR_PGIBUG_B)
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
    call atlas__output__write_fieldset(this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A, &
      config%CPTR_PGIBUG_B)
  else
    opt_config = atlas_Config()
    call atlas__output__write_fieldset(this%CPTR_PGIBUG_A,fieldset%CPTR_PGIBUG_A, &
      opt_config%CPTR_PGIBUG_B)
    call opt_config%final()
  endif
end subroutine

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Output__final_auto(this)
  type(atlas_Output), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_Output__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

end module atlas_output_module
