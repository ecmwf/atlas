
module atlas_MeshGenerator_module


use, intrinsic :: iso_c_binding, only: c_ptr
use atlas_c_interop, only: c_str
use atlas_refcounted_module, only: atlas_refcounted
use atlas_Grid_module, only: atlas_Grid
use atlas_GridDistribution_module, only: atlas_GridDistribution
use atlas_Mesh_module, only: atlas_Mesh
use atlas_Config_module, only: atlas_Config

implicit none

private :: c_ptr
private :: c_str
private :: atlas_refcounted
private :: atlas_Grid
private :: atlas_GridDistribution
private :: atlas_Mesh
private :: atlas_Config

public :: atlas_MeshGenerator
public :: atlas_meshgenerator_Structured

private

!-----------------------------
! atlas_Mesh                 !
!-----------------------------

!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_MeshGenerator

! Purpose :
! -------
!   *Nabla* :

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_MeshGenerator__delete
  procedure, public :: copy => atlas_MeshGenerator__copy
  procedure, public :: generate => atlas_MeshGenerator__generate

END TYPE atlas_MeshGenerator

interface atlas_MeshGenerator
  module procedure atlas_MeshGenerator__cptr
  module procedure atlas_MeshGenerator__name_config
end interface

interface atlas_meshgenerator_Structured
  module procedure atlas_meshgenerator_Structured__config
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


function atlas_MeshGenerator__cptr(cptr) result(MeshGenerator)
  type(atlas_MeshGenerator) :: MeshGenerator
  type(c_ptr), intent(in) :: cptr
  call MeshGenerator%reset_c_ptr( cptr )
end function

function atlas_MeshGenerator__name_config(name,config) result(MeshGenerator)
  use atlas_MeshGenerator_c_binding
  type(atlas_MeshGenerator) :: meshgenerator
  character(len=*), intent(in) :: name
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    meshgenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str(name),config%c_ptr()))
  else
    opt_config = atlas_Config()
    meshgenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str(name),opt_config%c_ptr()))
    call opt_config%final()
  endif
  call meshgenerator%return()
end function

function atlas_meshgenerator_Structured__config(config) result(meshgenerator)
  use atlas_MeshGenerator_c_binding
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Config), intent(in), optional :: config
  type(atlas_Config) :: opt_config
  if( present(config) ) then
    meshgenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str("Structured"),config%c_ptr()))
  else
    opt_config = atlas_Config()
    meshgenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str("Structured"),opt_config%c_ptr()))
    call opt_config%final()
  endif
  call meshgenerator%return()
end function


subroutine atlas_MeshGenerator__delete(this)
  use atlas_MeshGenerator_c_binding
  class(atlas_MeshGenerator), intent(inout) :: this
  if ( .not. this%is_null() ) then
    call atlas__MeshGenerator__delete(this%c_ptr())
  endif
  call this%reset_c_ptr()
end subroutine atlas_MeshGenerator__delete


subroutine atlas_MeshGenerator__copy(this,obj_in)
  class(atlas_MeshGenerator), intent(inout) :: this
  class(atlas_RefCounted), target, intent(in) :: obj_in
end subroutine

function atlas_MeshGenerator__generate(this,grid,distribution) result(mesh)
   use atlas_MeshGenerator_c_binding
   type(atlas_Mesh) :: mesh
   class(atlas_MeshGenerator), intent(in) :: this
   class(atlas_Grid), intent(in) :: grid
   class(atlas_GridDistribution), intent(in), optional :: distribution

   if( present(distribution) ) then
     mesh = atlas_Mesh( atlas__MeshGenerator__generate__grid_griddist(this%c_ptr(),grid%c_ptr(),distribution%c_ptr()) )
   else
     mesh = atlas_Mesh( atlas__MeshGenerator__generate__grid(this%c_ptr(),grid%c_ptr()) )
   endif
   call mesh%return()
end function


! -----------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------

end module atlas_MeshGenerator_module
