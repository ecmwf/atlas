
module atlas_MeshGenerator_module


use fckit_refcounted_module, only: fckit_refcounted

implicit none

private :: fckit_refcounted

public :: atlas_MeshGenerator

private

!-----------------------------!
! atlas_MeshGenerator         !
!-----------------------------!

!------------------------------------------------------------------------------
TYPE, extends(fckit_refcounted) :: atlas_MeshGenerator
contains
  procedure, public :: delete => atlas_MeshGenerator__delete
  procedure, public :: copy => atlas_MeshGenerator__copy
  procedure, public :: generate => atlas_MeshGenerator__generate
END TYPE atlas_MeshGenerator

interface atlas_MeshGenerator
  module procedure atlas_MeshGenerator__cptr
  module procedure atlas_MeshGenerator__config
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


function atlas_MeshGenerator__cptr(cptr) result(MeshGenerator)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_MeshGenerator) :: MeshGenerator
  type(c_ptr), intent(in) :: cptr
  call MeshGenerator%reset_c_ptr( cptr )
end function


function atlas_MeshGenerator__config(config) result(meshgenerator)
  use fckit_c_interop_module, only: c_str
  use atlas_MeshGenerator_c_binding
  use atlas_Config_module, only: atlas_Config
  type(atlas_MeshGenerator) :: meshgenerator
  type(atlas_Config), intent(in), optional :: config
  character(len=:), allocatable :: meshgenerator_type
  if( present(config) ) then
    if( .not. config%get("type",meshgenerator_type) ) then
       allocate(meshgenerator_type, source='structured')
    endif
    meshgenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create(c_str(meshgenerator_type),config%c_ptr()))
  else
    meshgenerator = atlas_MeshGenerator__cptr(atlas__MeshGenerator__create_noconfig('structured'))
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
  class(fckit_refcounted), target, intent(in) :: obj_in
end subroutine

function atlas_MeshGenerator__generate(this,grid,distribution) result(mesh)
   use atlas_MeshGenerator_c_binding
   use atlas_Grid_module, only: atlas_Grid
   use atlas_GridDistribution_module, only: atlas_GridDistribution
   use atlas_Mesh_module, only: atlas_Mesh
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
