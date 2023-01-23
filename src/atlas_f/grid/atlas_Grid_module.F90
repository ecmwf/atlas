! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_Grid_module

use fckit_owned_object_module, only: fckit_owned_object
use atlas_Config_module, only: atlas_Config
use atlas_Projection_module, only : atlas_Projection
use atlas_Domain_module, only : atlas_LonLatRectangularDomain
use atlas_kinds_module, only : ATLAS_KIND_IDX, ATLAS_KIND_GIDX
use, intrinsic :: iso_c_binding, only : c_ptr

implicit none

private :: fckit_owned_object
private :: atlas_Config
private :: atlas_Projection
private :: c_ptr

public :: atlas_Grid
public :: atlas_UnstructuredGrid
public :: atlas_StructuredGrid
public :: atlas_GaussianGrid
public :: atlas_ReducedGaussianGrid
public :: atlas_RegularGaussianGrid
public :: atlas_RegularLonLatGrid
public :: atlas_ShiftedLonLatGrid
public :: atlas_ShiftedLonGrid
public :: atlas_ShiftedLatGrid
public :: atlas_RegionalGrid

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_Grid

! Purpose :
! -------
!   *atlas_Grid* : Object Grid specifications for Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: size => atlas_Grid__size
  procedure :: spec => atlas_Grid__spec
  procedure :: uid
  procedure :: lonlat_bounding_box

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_Grid__final_auto
#endif

END TYPE atlas_Grid

interface atlas_Grid
  module procedure atlas_Grid__ctor_id
  module procedure atlas_Grid__ctor_config
  module procedure atlas_Grid__ctor_cptr
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_Grid) :: atlas_UnstructuredGrid

! Purpose :
! -------
!   *atlas_UnstructuredGrid* : Object Grid specifications for Reduced Grids

! Methods :
! -------

! Author :
! ------
!   8-Jun-2019 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_UnstructuredGrid__final_auto
#endif
END TYPE atlas_UnstructuredGrid

interface atlas_UnstructuredGrid
  module procedure atlas_UnstructuredGrid__ctor_points
  module procedure atlas_UnstructuredGrid__ctor_config
  module procedure atlas_UnstructuredGrid__ctor_cptr
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_Grid) :: atlas_StructuredGrid

! Purpose :
! -------
!   *atlas_StructuredGrid* : Object Grid specifications for Reduced Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: ny        => Structured__ny
  procedure, private   :: nx_int32 => Structured__nx_int32
  procedure, private   :: nx_int64 => Structured__nx_int64
  generic   :: nx => nx_int32, nx_int64
  procedure, private :: Structured__index_int32
  procedure, private :: Structured__index_int64
  generic   :: index => Structured__index_int32, Structured__index_int64
  procedure, private :: Structured__index2ij_int32
  procedure, private :: Structured__index2ij_int64
  generic   :: index2ij => Structured__index2ij_int32, Structured__index2ij_int64
  procedure, private :: Structured__ij_int32
  procedure, private :: Structured__ij_int64
  generic   :: ij => Structured__ij_int32, Structured__ij_int64
  procedure :: nx_array  => Structured__nx_array
  procedure :: nxmin     => Structured__nxmin
  procedure :: nxmax     => Structured__nxmax
  procedure :: y_array   => Structured__y_array
  procedure, private :: x_32  => Structured__x_32
  procedure, private :: x_64  => Structured__x_64
  generic  :: x         => x_32, x_64
  procedure, private :: y_32         => Structured__y_32
  procedure, private :: y_64         => Structured__y_64
  generic :: y         => y_32, y_64
  procedure, private :: xy_32        => Structured__xy_32
  procedure, private :: xy_64        => Structured__xy_64
  generic :: xy        => xy_32, xy_64
  procedure, private :: lonlat_32    => Structured__lonlat_32
  procedure, private :: lonlat_64    => Structured__lonlat_64
  generic :: lonlat    => lonlat_32, lonlat_64
  procedure :: reduced   => Structured__reduced
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_StructuredGrid__final_auto
#endif

END TYPE atlas_StructuredGrid

interface atlas_StructuredGrid
  module procedure atlas_StructuredGrid__ctor_id
  module procedure atlas_StructuredGrid__ctor_config
  module procedure atlas_StructuredGrid__ctor_cptr
end interface

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
TYPE, extends(atlas_StructuredGrid) :: atlas_GaussianGrid

! Purpose :
! -------
!   *atlas_ReducedGaussianGrid* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
procedure :: N         => Gaussian__N

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_GaussianGrid__final_auto
#endif

END TYPE atlas_GaussianGrid

interface atlas_GaussianGrid
  module procedure atlas_GaussianGrid__ctor_id
end interface

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_ReducedGaussianGrid

! Purpose :
! -------
!   *atlas_ReducedGaussianGrid* : Object Grid specifications for Reduced Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
procedure :: N         => ReducedGaussian__N

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_ReducedGaussianGrid__final_auto
#endif

END TYPE atlas_ReducedGaussianGrid

interface atlas_ReducedGaussianGrid
  module procedure atlas_ReducedGaussianGrid__ctor_int32
  module procedure atlas_ReducedGaussianGrid__ctor_int64
  module procedure atlas_ReducedGaussianGrid__ctor_projection_int32
  module procedure atlas_ReducedGaussianGrid__ctor_projection_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_RegularGaussianGrid

! Purpose :
! -------
!   *atlas_RegularGaussianGrid* : Object Grid specifications for Regular Gaussian Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
procedure :: N         => RegularGaussian__N

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_RegularGaussianGrid__final_auto
#endif

END TYPE atlas_RegularGaussianGrid

interface atlas_RegularGaussianGrid
  module procedure atlas_RegularGaussianGrid__ctor_int32
  module procedure atlas_RegularGaussianGrid__ctor_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_RegularLonLatGrid

! Purpose :
! -------
!   *atlas_RegularLonLatGrid* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_RegularLonLatGrid__final_auto
#endif
END TYPE atlas_RegularLonLatGrid

interface atlas_RegularLonLatGrid
  module procedure atlas_grid_RegularLonLat__ctor_int32
  module procedure atlas_grid_RegularLonLat__ctor_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_ShiftedLonLatGrid

! Purpose :
! -------
!   *atlas_ShiftedLonLatGrid* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_ShiftedLonLatGrid__final_auto
#endif
END TYPE atlas_ShiftedLonLatGrid

interface atlas_ShiftedLonLatGrid
  module procedure atlas_grid_ShiftedLonLat__ctor_int32
  module procedure atlas_grid_ShiftedLonLat__ctor_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_ShiftedLonGrid

! Purpose :
! -------
!   *atlas_ShiftedLonGrid* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_ShiftedLonGrid__final_auto
#endif
END TYPE atlas_ShiftedLonGrid

interface atlas_ShiftedLonGrid
  module procedure atlas_grid_ShiftedLon__ctor_int32
  module procedure atlas_grid_ShiftedLon__ctor_int64
end interface

!------------------------------------------------------------------------------

TYPE, extends(atlas_StructuredGrid) :: atlas_ShiftedLatGrid

! Purpose :
! -------
!   *atlas_ShiftedLatGrid* : Object Grid specifications for LonLat Grids

! Methods :
! -------

! Author :
! ------
!   9-Oct-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_ShiftedLatGrid__final_auto
#endif
END TYPE atlas_ShiftedLatGrid

interface atlas_ShiftedLatGrid
  module procedure atlas_grid_ShiftedLat__ctor_int32
  module procedure atlas_grid_ShiftedLat__ctor_int64
end interface

!------------------------------------------------------------------------------

interface atlas_RegionalGrid
    module procedure atlas_RegionalGrid_ctor_int32
    module procedure atlas_RegionalGrid_ctor_int64
    module procedure atlas_RegionalGrid_ctor_nwse_int32
    module procedure atlas_RegionalGrid_ctor_nwse_int64
    module procedure atlas_RegionalGrid_ctor_increments_int32
    module procedure atlas_RegionalGrid_ctor_increments_int64
end interface

!------------------------------------------------------------------------------
interface c_idx
  module procedure c_idx_32
  module procedure c_idx_64
end interface

interface c_gidx
  module procedure c_gidx_32
  module procedure c_gidx_64
end interface

!------------------------------------------------------------------------------

!========================================================
contains
!========================================================

pure function c_idx_32(f_idx) result(c_idx)
    use, intrinsic :: iso_c_binding, only : c_long
    integer(ATLAS_KIND_IDX) :: c_idx
    integer(c_long), intent(in) :: f_idx
    c_idx = int(f_idx,ATLAS_KIND_IDX) - 1_ATLAS_KIND_IDX
end function

pure function c_idx_64(f_idx) result(c_idx)
    use, intrinsic :: iso_c_binding, only : c_long, c_int
    integer(ATLAS_KIND_IDX) :: c_idx
    integer(c_int), intent(in) :: f_idx
    c_idx = int(f_idx,ATLAS_KIND_IDX) - 1_ATLAS_KIND_IDX
end function

pure function c_gidx_32(f_gidx) result(c_gidx)
    use, intrinsic :: iso_c_binding, only : c_long
    integer(ATLAS_KIND_GIDX) :: c_gidx
    integer(c_long), intent(in) :: f_gidx
    c_gidx = int(f_gidx,ATLAS_KIND_GIDX) - 1_ATLAS_KIND_GIDX
end function

pure function c_gidx_64(f_gidx) result(c_gidx)
    use, intrinsic :: iso_c_binding, only : c_int
    integer(ATLAS_KIND_GIDX) :: c_gidx
    integer(c_int), intent(in) :: f_gidx
    c_gidx = int(f_gidx,ATLAS_KIND_GIDX) - 1_ATLAS_KIND_GIDX
end function

! -----------------------------------------------------------------------------
! Destructor

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_Grid__final_auto(this)
  type(atlas_Grid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_UnstructuredGrid__final_auto(this)
  type(atlas_UnstructuredGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_StructuredGrid__final_auto(this)
  type(atlas_StructuredGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_GaussianGrid__final_auto(this)
  type(atlas_GaussianGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_ReducedGaussianGrid__final_auto(this)
  type(atlas_ReducedGaussianGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_RegularLonLatGrid__final_auto(this)
  type(atlas_RegularLonLatGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_ShiftedLonLatGrid__final_auto(this)
  type(atlas_ShiftedLonLatGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_ShiftedLonGrid__final_auto(this)
  type(atlas_ShiftedLonGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_ShiftedLatGrid__final_auto(this)
  type(atlas_ShiftedLatGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

ATLAS_FINAL subroutine atlas_RegularGaussianGrid__final_auto(this)
  type(atlas_RegularGaussianGrid), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif


! -----------------------------------------------------------------------------
! Constructors

function atlas_Grid__ctor_id(identifier) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_Grid) :: this
  character(len=*), intent(in) :: identifier
  call this%reset_c_ptr( atlas__grid__Structured(c_str(identifier)) )
  call this%return()
end function

function atlas_Grid__ctor_config(config) result(this)
  use atlas_grid_Structured_c_binding
  type(atlas_Grid) :: this
  type(atlas_Config), intent(in) :: config
  call this%reset_c_ptr( atlas__grid__Structured__config(config%CPTR_PGIBUG_B) )
  call this%return()
end function

function atlas_Grid__ctor_cptr(cptr) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_Grid) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

! -----------------------------------------------------------------------------

function atlas_UnstructuredGrid__ctor_config(config) result(this)
  use atlas_grid_unstructured_c_binding
  type(atlas_UnstructuredGrid) :: this
  type(atlas_Config), intent(in) :: config
  call this%reset_c_ptr( atlas__grid__Unstructured__config(config%CPTR_PGIBUG_B) )
  call this%return()
end function

function atlas_UnstructuredGrid__ctor_cptr(cptr) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_unstructured_c_binding
  type(atlas_UnstructuredGrid) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_UnstructuredGrid__ctor_points( xy ) result(this)
  use, intrinsic :: iso_c_binding, only : c_double, c_int
  use fckit_c_interop_module, only : c_str
  use fckit_array_module, only : array_strides, array_view1d
  use atlas_grid_unstructured_c_binding
  type(atlas_UnstructuredGrid) :: this
  real(c_double), intent(in) :: xy(:,:)
  integer(c_int) :: shapef(2)
  integer(c_int) :: stridesf(2)
  real(c_double), pointer :: xy1d(:)
  xy1d => array_view1d(xy)
  shapef = shape(xy)
  stridesf = array_strides(xy)
  call this%reset_c_ptr( atlas__grid__Unstructured__points( xy1d, shapef, stridesf ) )
  call this%return()
end function

! -----------------------------------------------------------------------------

function atlas_StructuredGrid__ctor_id(identifier) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_StructuredGrid) :: this
  character(len=*), intent(in) :: identifier
  call this%reset_c_ptr( atlas__grid__Structured(c_str(identifier)) )
  call this%return()
end function

function atlas_StructuredGrid__ctor_config(config) result(this)
  use atlas_grid_Structured_c_binding
  type(atlas_StructuredGrid) :: this
  type(atlas_Config), intent(in) :: config
  call this%reset_c_ptr( atlas__grid__Structured__config(config%CPTR_PGIBUG_B) )
  call this%return()
end function

function atlas_StructuredGrid__ctor_cptr(cptr) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_StructuredGrid) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

!-----------------------------------------------------------------------------

function atlas_GaussianGrid__ctor_id(identifier) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_grid_Structured_c_binding
  type(atlas_GaussianGrid) :: this
  character(len=*), intent(in) :: identifier
  call this%reset_c_ptr( atlas__grid__Structured(c_str(identifier)) )
  call this%return()
end function

! -----------------------------------------------------------------------------

function atlas_RegularGaussianGrid__ctor_int32(N) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularGaussianGrid) :: this
  integer(c_int), intent(in) :: N
  call this%reset_c_ptr( atlas__grid__regular__RegularGaussian(int(N,c_long)) )
  call this%return()
end function

function atlas_RegularGaussianGrid__ctor_int64(N) result(this)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularGaussianGrid) :: this
  integer(c_long), intent(in) :: N
  call this%reset_c_ptr( atlas__grid__regular__RegularGaussian(int(N,c_long)) )
  call this%return()
end function

!-----------------------------------------------------------------------------

function atlas_ReducedGaussianGrid__ctor_int32(nx) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: this
  integer(c_int), intent(in)  :: nx(:)
  call this%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_int( nx, int(size(nx),c_long) ) )
   call this%return()
end function

function atlas_ReducedGaussianGrid__ctor_int64(nx) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: this
  integer(c_long), intent(in)  :: nx(:)
  call this%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_long( nx, int(size(nx),c_long) ) )
  call this%return()
end function

function atlas_ReducedGaussianGrid__ctor_projection_int32(nx, projection) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: this
  integer(c_int), intent(in)  :: nx(:)
  type(atlas_Projection), intent(in) :: projection
  call this%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_int_projection( nx, int(size(nx),c_long), projection%CPTR_PGIBUG_A ) )
   call this%return()
end function

function atlas_ReducedGaussianGrid__ctor_projection_int64(nx, projection) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ReducedGaussianGrid) :: this
  integer(c_long), intent(in)  :: nx(:)
  type(atlas_Projection), intent(in) :: projection
  call this%reset_c_ptr( &
    & atlas__grid__reduced__ReducedGaussian_long_projection( nx, int(size(nx),c_long), projection%CPTR_PGIBUG_A ) )
  call this%return()
end function

!-----------------------------------------------------------------------------

function atlas_grid_RegularLonLat__ctor_int32(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularLonLatGrid) :: this
  integer(c_int), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__RegularLonLat(int(nlon,c_long),int(nlat,c_long)) )
  call this%return()
end function

function atlas_grid_RegularLonLat__ctor_int64(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_RegularLonLatGrid) :: this
  integer(c_long), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__RegularLonLat( nlon, nlat ) )
  call this%return()
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLonLat__ctor_int32(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ShiftedLonLatGrid) :: this
  integer(c_int), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__ShiftedLonLat(int(nlon,c_long),int(nlat,c_long)) )
  call this%return()
end function

function atlas_grid_ShiftedLonLat__ctor_int64(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ShiftedLonLatGrid) :: this
  integer(c_long), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__ShiftedLonLat( nlon, nlat ) )
  call this%return()
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLon__ctor_int32(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ShiftedLonGrid) :: this
  integer(c_int), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__ShiftedLon(int(nlon,c_long),int(nlat,c_long)) )
  call this%return()
end function

function atlas_grid_ShiftedLon__ctor_int64(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ShiftedLonGrid) :: this
  integer(c_long), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__ShiftedLon( nlon, nlat ) )
  call this%return()
end function

!-----------------------------------------------------------------------------

function atlas_grid_ShiftedLat__ctor_int32(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ShiftedLatGrid) :: this
  integer(c_int), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__ShiftedLat(int(nlon,c_long),int(nlat,c_long)) )
  call this%return()
end function

function atlas_grid_ShiftedLat__ctor_int64(nlon,nlat) result(this)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  type(atlas_ShiftedLatGrid) :: this
  integer(c_long), intent(in) :: nlon, nlat
  call this%reset_c_ptr( atlas__grid__regular__ShiftedLat( nlon, nlat ) )
  call this%return()
end function

! -----------------------------------------------------------------------------
! Structured members

function atlas_Grid__size(this) result(npts)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_grid_Grid_c_binding
  class(atlas_Grid), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: npts
  npts = int(atlas__grid__Grid__size(this%CPTR_PGIBUG_A),ATLAS_KIND_IDX)
end function

function atlas_Grid__spec(this) result(spec)
  use atlas_grid_Grid_c_binding
  class(atlas_Grid), intent(in) :: this
  type(atlas_Config) :: spec
  spec = atlas_Config( atlas__grid__Grid__spec(this%CPTR_PGIBUG_A) )
  call spec%return ()
end function

function uid(this)
  use atlas_grid_Grid_c_binding
  use fckit_c_interop_module, only : c_ptr_to_string, c_ptr_free
  use, intrinsic :: iso_c_binding, only : c_ptr
  class(atlas_Grid), intent(in) :: this
  character(len=:), allocatable :: uid
  type(c_ptr) :: uid_c_str
  integer :: size
  call atlas__grid__Grid__uid(this%CPTR_PGIBUG_A, uid_c_str, size )
  uid = c_ptr_to_string(uid_c_str)
  call c_ptr_free(uid_c_str)
end function

function lonlat_bounding_box (this) result (bb)
  use, intrinsic :: iso_c_binding, only : c_double
  use atlas_grid_Grid_c_binding
  class(atlas_Grid), intent(in) :: this
  type(atlas_LonLatRectangularDomain) :: bb
  bb = atlas_LonLatRectangularDomain(atlas__grid__Grid__lonlat_bounding_box (this%CPTR_PGIBUG_A))
  call bb%return()
end function


function Gaussian__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_GaussianGrid), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Gaussian__N(this%CPTR_PGIBUG_A)
end function

function ReducedGaussian__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_ReducedGaussianGrid), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Gaussian__N(this%CPTR_PGIBUG_A)
end function

function RegularGaussian__N(this) result(N)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_RegularGaussianGrid), intent(in) :: this
  integer(c_long) :: N
  N = atlas__grid__Gaussian__N(this%CPTR_PGIBUG_A)
end function

function Structured__ny(this) result(ny)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: ny
  ny = atlas__grid__Structured__ny(this%CPTR_PGIBUG_A)
end function

subroutine Structured__index2ij_int32(this, gidx, i, j) 
  use, intrinsic :: iso_c_binding, only: c_int32_t
  use atlas_grid_Structured_c_binding
  integer(c_int32_t), intent (in) :: gidx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(ATLAS_KIND_IDX), intent(out) :: i, j
  call atlas__grid__Structured__index2ij(this%CPTR_PGIBUG_A, c_gidx(gidx), i, j)
  i = i + 1
  j = j + 1
end subroutine

subroutine Structured__index2ij_int64(this, gidx, i, j) 
  use, intrinsic :: iso_c_binding, only: c_long, c_int64_t
  use atlas_grid_Structured_c_binding
  integer(c_int64_t), intent (in) :: gidx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(ATLAS_KIND_IDX), intent(out) :: i, j
  call atlas__grid__Structured__index2ij(this%CPTR_PGIBUG_A, c_gidx(gidx), i, j)
  i = i + 1
  j = j + 1
end subroutine

function Structured__ij_int32(this, gidx) result(ij)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_grid_Structured_c_binding
  integer(c_int), intent (in) :: gidx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: ij (2)
  call atlas__grid__Structured__index2ij(this%CPTR_PGIBUG_A, c_gidx(gidx), ij(1), ij(2))
  ij = ij + 1
end function

function Structured__ij_int64(this, gidx) result(ij)
  use, intrinsic :: iso_c_binding, only: c_int, c_long
  use atlas_grid_Structured_c_binding
  integer(c_long), intent (in) :: gidx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: ij (2)
  call atlas__grid__Structured__index2ij(this%CPTR_PGIBUG_A, c_gidx(gidx), ij(1), ij(2))
  ij = ij + 1
end function

function Structured__index_int32(this, i, j) result(gidx)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_grid_Structured_c_binding
  integer(ATLAS_KIND_IDX) :: gidx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int), intent(in) :: i, j
  gidx = int(1 + atlas__grid__Structured__index(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j) ), ATLAS_KIND_IDX )
end function

function Structured__index_int64(this, i, j) result(gidx)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  integer(ATLAS_KIND_IDX) :: gidx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long), intent(in) :: i, j
  gidx = int( 1 + atlas__grid__Structured__index(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j) ), ATLAS_KIND_IDX )
end function

function Structured__nx_int32(this, j) result(nx)
  use, intrinsic :: iso_c_binding, only: c_int
  use atlas_grid_Structured_c_binding
  integer(ATLAS_KIND_IDX) :: nx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int), intent(in) :: j
  nx = atlas__grid__Structured__nx(this%CPTR_PGIBUG_A, c_idx(j) )
end function

function Structured__nx_int64(this, j) result(nx)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  integer(ATLAS_KIND_IDX) :: nx
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long), intent(in) :: j
  nx = atlas__grid__Structured__nx(this%CPTR_PGIBUG_A, c_idx(j) )
end function

function Structured__reduced(this) result(reduced)
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in) :: this
  logical :: reduced
  if( atlas__grid__Structured__reduced(this%CPTR_PGIBUG_A) == 1 ) then
    reduced = .true.
  else
    reduced = .false.
  endif
end function

function Structured__nx_array(this) result(nx)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_long, c_f_pointer
  class(atlas_StructuredGrid), intent(in) :: this
  integer(ATLAS_KIND_IDX), pointer        :: nx(:)
  type   (c_ptr)                          :: nx_c_ptr
  integer(ATLAS_KIND_IDX)                 :: nx_size
  call atlas__grid__Structured__nx_array(this%CPTR_PGIBUG_A, nx_c_ptr, nx_size)
  call c_f_pointer(nx_c_ptr , nx , (/nx_size/))
end function

function Structured__nxmax(this) result(nxmax)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  integer(c_long)                          :: nxmax
  nxmax = atlas__grid__Structured__nxmax(this%CPTR_PGIBUG_A)
end function

function Structured__nxmin(this) result(nxmin)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  integer(c_long)                          :: nxmin
  nxmin = atlas__grid__Structured__nxmin(this%CPTR_PGIBUG_A)
end function

function Structured__y_32(this, j) result(y)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  real(c_double) :: y
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int),              intent(in) :: j
  y = atlas__grid__Structured__y(this%CPTR_PGIBUG_A, c_idx(j) )
end function

function Structured__y_64(this, j) result(y)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: y
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long),             intent(in) :: j
  y = atlas__grid__Structured__y(this%CPTR_PGIBUG_A, c_idx(j) )
end function

function Structured__y_array(this) result(y)
  use atlas_grid_Structured_c_binding
  use, intrinsic :: iso_c_binding , only : c_double, c_f_pointer
  class(atlas_StructuredGrid), intent(in) :: this
  real   (c_double)       , pointer    :: y(:)
  type   (c_ptr)                       :: y_c_ptr
  integer(ATLAS_KIND_IDX)              :: y_size
  call atlas__grid__Structured__y_array(this%CPTR_PGIBUG_A, &
      & y_c_ptr, y_size)
  call c_f_pointer (y_c_ptr, y, (/y_size/))
end function

function Structured__x_32(this, i,j) result(x)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  real(c_double) :: x
  integer(c_int) :: i,j
  x = atlas__grid__Structured__x(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j))
end function

function Structured__x_64(this, i,j) result(x)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  class(atlas_StructuredGrid), intent(in)  :: this
  real(c_double) :: x
  integer(c_long) :: i,j
  x = atlas__grid__Structured__x(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j))
end function

function Structured__xy_32(this, i,j) result(xy)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  real(c_double) :: xy(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int) , intent(in) :: i,j
  call atlas__grid__Structured__xy(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j), xy)
end function

function Structured__xy_64(this, i,j) result(xy)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: xy(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) , intent(in) :: i,j
  call atlas__grid__Structured__xy(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j), xy)
end function

function Structured__lonlat_32(this, i,j) result(lonlat)
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use atlas_grid_Structured_c_binding
  real(c_double) :: lonlat(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_int) , intent(in) :: i,j
  call atlas__grid__Structured__lonlat(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j), lonlat)
end function

function Structured__lonlat_64(this, i,j) result(lonlat)
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use atlas_grid_Structured_c_binding
  real(c_double) :: lonlat(2)
  class(atlas_StructuredGrid), intent(in) :: this
  integer(c_long) , intent(in) :: i,j
  call atlas__grid__Structured__lonlat(this%CPTR_PGIBUG_A, c_idx(i), c_idx(j), lonlat)
end function

! ----------------------------------------------------------------------------------------

function atlas_RegionalGrid_ctor_int32( nx, ny, xy_min, xy_max, projection, y_numbering ) result(this)
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    type(atlas_StructuredGrid) :: this
    integer(c_int), intent(in) :: nx, ny
    real(c_double), intent(in) :: xy_min(2), xy_max(2)
    type(atlas_Projection), intent(in), optional :: projection
    integer(c_int), intent(in), optional :: y_numbering
    type(atlas_Config) :: config

    config = atlas_Config()
    call config%set("type","regional")
    call config%set("nx",nx)
    call config%set("ny",ny)
    call config%set("xmin",xy_min(1))
    call config%set("ymin",xy_min(2))
    call config%set("xmax",xy_max(1))
    call config%set("ymax",xy_max(2))
    if( present(projection) ) then
      call config%set("projection",projection%spec())
    endif
    if( present(y_numbering) ) then
      call config%set("y_numbering",y_numbering)
    endif
    this = atlas_StructuredGrid(config)
    call config%final()
    call this%return()
end function

function atlas_RegionalGrid_ctor_int64( nx, ny, xy_min, xy_max, projection, y_numbering ) result(this)
    use, intrinsic :: iso_c_binding, only : c_int, c_long, c_double
    type(atlas_StructuredGrid) :: this
    integer(c_long), intent(in) :: nx, ny
    real(c_double), intent(in) :: xy_min(2), xy_max(2)
    type(atlas_Projection), intent(in), optional :: projection
    integer(c_int), intent(in), optional :: y_numbering
    type(atlas_Config) :: config

    config = atlas_Config()
    call config%set("type","regional")
    call config%set("nx",nx)
    call config%set("ny",ny)
    call config%set("xmin",xy_min(1))
    call config%set("ymin",xy_min(2))
    call config%set("xmax",xy_max(1))
    call config%set("ymax",xy_max(2))
    if( present(projection) ) then
      call config%set("projection",projection%spec())
    endif
    if( present(y_numbering) ) then
      call config%set("y_numbering",y_numbering)
    endif
    this = atlas_StructuredGrid(config)
    call config%final()
    call this%return()
end function


function atlas_RegionalGrid_ctor_nwse_int32( nx, ny, north, west, south, east, projection, y_numbering ) result(this)
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    type(atlas_StructuredGrid) :: this
    integer(c_int), intent(in) :: nx, ny
    real(c_double), intent(in) :: north, west, south, east
    type(atlas_Projection), intent(in), optional :: projection
    integer(c_int), intent(in), optional :: y_numbering
    type(atlas_Config) :: config

    config = atlas_Config()
    call config%set("type","regional")
    call config%set("nx",nx)
    call config%set("ny",ny)
    call config%set("north",north)
    call config%set("west",west)
    call config%set("south",south)
    call config%set("east",east)
    if( present(projection) ) then
      call config%set("projection",projection%spec())
    endif
    if( present(y_numbering) ) then
      call config%set("y_numbering",y_numbering)
    endif
    this = atlas_StructuredGrid(config)
    call config%final()
    call this%return()
end function

function atlas_RegionalGrid_ctor_nwse_int64( nx, ny, north, west, south, east, projection, y_numbering ) result(this)
    use, intrinsic :: iso_c_binding, only : c_int, c_long, c_double
    type(atlas_Grid) :: this
    integer(c_long), intent(in) :: nx, ny
    real(c_double), intent(in) :: north, west, south, east
    type(atlas_Projection), intent(in), optional :: projection
    integer(c_int), intent(in), optional :: y_numbering
    type(atlas_Config) :: config

    config = atlas_Config()
    call config%set("type","regional")
    call config%set("nx",nx)
    call config%set("ny",ny)
    call config%set("north",north)
    call config%set("west",west)
    call config%set("south",south)
    call config%set("east",east)
    if( present(projection) ) then
      call config%set("projection",projection%spec())
    endif
    if( present(y_numbering) ) then
      call config%set("y_numbering",y_numbering)
    endif
    this = atlas_StructuredGrid(config)
    call config%final()
    call this%return()
end function


function atlas_RegionalGrid_ctor_increments_int32( nx, ny, dx, dy, xy_min, projection, y_numbering ) result(this)
    use iso_c_binding, only : c_int, c_double
    type(atlas_StructuredGrid) :: this
    integer(c_int), intent(in) :: nx
    integer(c_int), intent(in) :: ny
    real(c_double), intent(in) :: xy_min(2)
    real(c_double), intent(in) :: dx
    real(c_double), intent(in) :: dy
    type(atlas_Projection), intent(in), optional :: projection
    integer(c_int), intent(in), optional :: y_numbering

    type(atlas_Config) :: config

    config = atlas_Config()
    call config%set("type","regional")
    call config%set("nx",nx)
    call config%set("ny",ny)
    call config%set("xmin",xy_min(1))
    call config%set("ymin",xy_min(2))
    call config%set("dx",dx)
    call config%set("dy",dy)
    if( present(projection) ) then
      call config%set("projection",projection%spec())
    endif
    if( present(y_numbering) ) then
      call config%set("y_numbering",y_numbering)
    endif

    this = atlas_StructuredGrid(config)
    call config%final()
    call this%return()
end function

function atlas_RegionalGrid_ctor_increments_int64( nx, ny, dx, dy, xy_min, projection, y_numbering ) result(this)
    use iso_c_binding, only : c_int, c_long, c_double
    type(atlas_StructuredGrid) :: this
    integer(c_long), intent(in) :: nx
    integer(c_long), intent(in) :: ny
    real(c_double), intent(in) :: xy_min(2)
    real(c_double), intent(in) :: dx
    real(c_double), intent(in) :: dy
    type(atlas_Projection), intent(in), optional :: projection
    integer(c_int), intent(in), optional :: y_numbering

    type(atlas_Config) :: config

    config = atlas_Config()
    call config%set("type","regional")
    call config%set("nx",nx)
    call config%set("ny",ny)
    call config%set("xmin",xy_min(1))
    call config%set("ymin",xy_min(2))
    call config%set("dx",dx)
    call config%set("dy",dy)
    if( present(projection) ) then
      call config%set("projection",projection%spec())
    endif
    if( present(y_numbering) ) then
      call config%set("y_numbering",y_numbering)
    endif

    this = atlas_StructuredGrid(config)
    call config%final()
    call this%return()
end function

! ----------------------------------------------------------------------------------------

end module atlas_Grid_module
