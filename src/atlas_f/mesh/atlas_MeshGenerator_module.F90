! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_MeshGenerator_module


use fckit_owned_object_module, only: fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_MeshGenerator

private

!-----------------------------!
! atlas_MeshGenerator         !
!-----------------------------!

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_MeshGenerator
contains
  procedure, public :: generate => atlas_MeshGenerator__generate
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_MeshGenerator__final_auto
#endif

END TYPE atlas_MeshGenerator

interface atlas_MeshGenerator
  module procedure atlas_MeshGenerator__cptr
  module procedure atlas_MeshGenerator__config
  module procedure atlas_MeshGenerator__type
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


function atlas_MeshGenerator__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_MeshGenerator) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_MeshGenerator__type(type) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_MeshGenerator_c_binding
  type(atlas_MeshGenerator) :: this
  character(len=*), intent(in) :: type
  call this%reset_c_ptr( atlas__MeshGenerator__create_noconfig(c_str(type)) )
  call this%return()
end function

function atlas_MeshGenerator__config(config) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_MeshGenerator_c_binding
  use atlas_Config_module, only: atlas_Config
  type(atlas_MeshGenerator) :: this
  type(atlas_Config), intent(in), optional :: config
  character(len=:), allocatable :: meshgenerator_type
  if( present(config) ) then
    if( .not. config%get("type",meshgenerator_type) ) then
       meshgenerator_type='structured'
    endif
    call this%reset_c_ptr( atlas__MeshGenerator__create( &
      c_str(meshgenerator_type),config%CPTR_PGIBUG_B) )
  else
    call this%reset_c_ptr( atlas__MeshGenerator__create_noconfig(c_str('structured')) )
  endif
  call this%return()
end function

function atlas_MeshGenerator__generate(this,grid,distribution) result(mesh)
   use atlas_MeshGenerator_c_binding
   use atlas_Grid_module, only: atlas_Grid
   use atlas_GridDistribution_module, only: atlas_GridDistribution
   use atlas_Mesh_module, only: atlas_Mesh
   type(atlas_Mesh) :: mesh
   class(atlas_MeshGenerator), intent(in) :: this
   class(atlas_Grid), intent(in) :: grid
   class(atlas_GridDistribution), intent(in), optional :: distribution
   call mesh%reset_c_ptr() ! Somehow needed with PGI/16.7 and build-type "bit"
   if( present(distribution) ) then
     mesh = atlas_Mesh( atlas__MeshGenerator__generate__grid_griddist( &
       this%CPTR_PGIBUG_A,grid%CPTR_PGIBUG_A,distribution%CPTR_PGIBUG_A) )
   else
     mesh = atlas_Mesh( atlas__MeshGenerator__generate__grid( &
       this%CPTR_PGIBUG_A,grid%CPTR_PGIBUG_A) )
   endif
   call mesh%return()
end function

!-------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_MeshGenerator__final_auto(this)
  type(atlas_MeshGenerator), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_MeshGenerator__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_MeshGenerator_module
