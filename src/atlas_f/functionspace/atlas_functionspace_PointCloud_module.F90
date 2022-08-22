! (C) Copyright 2020 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_functionspace_PointCloud_module

use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
use atlas_functionspace_module, only : atlas_FunctionSpace
use atlas_Grid_module, only: atlas_Grid
use atlas_Field_module, only: atlas_Field
use atlas_kinds_module, only: ATLAS_KIND_GIDX

implicit none

private :: c_str, c_ptr_to_string, c_ptr_free
private :: atlas_FunctionSpace
private :: atlas_Field
private :: atlas_Grid
private :: ATLAS_KIND_GIDX

public :: atlas_functionspace_PointCloud

private

!------------------------------------------------------------------------------
TYPE, extends(atlas_FunctionSpace) :: atlas_functionspace_PointCloud

! Purpose :
! -------
!   *atlas_functionspace_PointCloud* : Interpretes point cloud fields

! Methods :
! -------

! Author :
! ------
!   March-2020 Benjamin Menetrier     *IRIT-JCSDA*

!------------------------------------------------------------------------------
contains

  procedure, public :: size => size_
  procedure, public :: lonlat

END TYPE atlas_functionspace_PointCloud

interface atlas_functionspace_PointCloud
  module procedure ctor_cptr
  module procedure ctor_lonlat
  module procedure ctor_lonlat_ghost
  module procedure ctor_grid
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

!------------------------------------------------------------------------------

function ctor_cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(atlas_functionspace_PointCloud) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

!------------------------------------------------------------------------------

function ctor_lonlat(lonlat) result(this)
  use atlas_functionspace_PointCloud_c_binding
  type(atlas_functionspace_PointCloud) :: this
  class(atlas_Field), intent(in) :: lonlat
  call this%reset_c_ptr( atlas__functionspace__PointCloud__new__lonlat( lonlat%CPTR_PGIBUG_A ) )
  call this%return()
end function

!------------------------------------------------------------------------------

function ctor_lonlat_ghost(lonlat,ghost) result(this)
  use atlas_functionspace_PointCloud_c_binding
  type(atlas_functionspace_PointCloud) :: this
  class(atlas_Field), intent(in) :: lonlat
  class(atlas_Field), intent(in) :: ghost
  call this%reset_c_ptr( atlas__functionspace__PointCloud__new__lonlat_ghost( lonlat%CPTR_PGIBUG_A, ghost%CPTR_PGIBUG_A  ) )
  call this%return()
end function

!------------------------------------------------------------------------------

function ctor_grid(grid) result(this)
  use atlas_functionspace_PointCloud_c_binding
  type(atlas_functionspace_PointCloud) :: this
  class(atlas_Grid), intent(in) :: grid
  call this%reset_c_ptr( atlas__functionspace__PointCloud__new__grid( grid%CPTR_PGIBUG_A ) )
  call this%return()
end function

!------------------------------------------------------------------------------

function size_(this)
  use atlas_functionspace_PointCloud_c_binding
  integer :: size_
  class(atlas_functionspace_PointCloud), intent(in) :: this
  size_ = atlas__fs__PointCloud__size(this%CPTR_PGIBUG_A)
end function

!------------------------------------------------------------------------------

function lonlat(this) result(field)
  use atlas_functionspace_PointCloud_c_binding
  type(atlas_Field) :: field
  class(atlas_functionspace_PointCloud), intent(in) :: this
  field = atlas_Field( atlas__fs__PointCloud__lonlat(this%CPTR_PGIBUG_A) )
  call field%return()
end function

!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_functionspace_PointCloud__final_auto(this)
  type(atlas_functionspace_PointCloud), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_functionspace_PointCloud__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

!------------------------------------------------------------------------------

end module atlas_functionspace_PointCloud_module
