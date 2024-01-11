! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_StencilComputer_module

use fckit_shared_object_module, only: fckit_shared_object, fckit_c_deleter
use atlas_kinds_module, only : ATLAS_KIND_IDX
use, intrinsic :: iso_c_binding, only : c_ptr, c_int

implicit none

private :: fckit_shared_object, fckit_c_deleter
private :: c_ptr, c_int
private :: ATLAS_KIND_IDX

public :: atlas_StructuredGrid_ComputeNorth
public :: atlas_StructuredGrid_ComputeWest
public :: atlas_StructuredGrid_ComputeStencil
public :: atlas_StructuredGrid_Stencil

private

!------------------------------------------------------------------------------
TYPE, extends(fckit_shared_object) :: atlas_StructuredGrid_ComputeNorth

! Purpose :
! -------
!   *atlas_StructuredGrid_ComputeNorth* : To compute latitude index north of latitude

! Methods :
! -------

! Author :
! ------
!   9-Jan-2024 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: atlas_StructuredGrid_ComputeNorth__execute_real32
  procedure, private :: atlas_StructuredGrid_ComputeNorth__execute_real64
  generic :: execute => atlas_StructuredGrid_ComputeNorth__execute_real32, &
                      & atlas_StructuredGrid_ComputeNorth__execute_real64

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_StructuredGrid_ComputeNorth__final_auto
#endif

END TYPE atlas_StructuredGrid_ComputeNorth

interface atlas_StructuredGrid_ComputeNorth
  module procedure atlas_StructuredGrid_ComputeNorth__ctor
end interface

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
TYPE, extends(fckit_shared_object) :: atlas_StructuredGrid_ComputeWest

! Purpose :
! -------
!   *atlas_StructuredGrid_ComputeWest* : To compute longitude index west of longitude at given latitude index

! Methods :
! -------

! Author :
! ------
!   9-Jan-2024 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

  procedure, private :: atlas_StructuredGrid_ComputeWest__execute_real32
  procedure, private :: atlas_StructuredGrid_ComputeWest__execute_real64
  generic :: execute => atlas_StructuredGrid_ComputeWest__execute_real32, &
                      & atlas_StructuredGrid_ComputeWest__execute_real64

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_StructuredGrid_ComputeWest__final_auto
#endif

END TYPE atlas_StructuredGrid_ComputeWest

interface atlas_StructuredGrid_ComputeWest
  module procedure atlas_StructuredGrid_ComputeWest__ctor
end interface

!------------------------------------------------------------------------------

integer(c_int), parameter, private :: STENCIL_MAX_WIDTH = 6 !-> maximum 6x6 stencil

TYPE atlas_StructuredGrid_Stencil
  integer(ATLAS_KIND_IDX) :: width
  integer(ATLAS_KIND_IDX) :: j_begin
  integer(ATLAS_KIND_IDX) :: i_begin(STENCIL_MAX_WIDTH) ! on stack
contains
  procedure, pass :: write => atlas_StructuredGrid_Stencil__write
  generic, public :: write(FORMATTED) => write
  procedure, public :: i => atlas_StructuredGrid_Stencil__i
  procedure, public :: j => atlas_StructuredGrid_Stencil__j
END TYPE atlas_StructuredGrid_Stencil

TYPE :: atlas_StructuredGrid_ComputeStencil
  integer(c_int) :: halo
  integer(ATLAS_KIND_IDX) :: stencil_width
  integer(ATLAS_KIND_IDX) :: stencil_offset
  type(atlas_StructuredGrid_ComputeNorth) :: compute_north
  type(atlas_StructuredGrid_ComputeWest)  :: compute_west
contains
  procedure :: setup => atlas_StructuredGrid_ComputeStencil__setup
  procedure :: execute => atlas_StructuredGrid_ComputeStencil__execute_real64

  procedure :: assignment_operator => atlas_StructuredGrid_ComputeStencil__assignment
  generic, public :: assignment(=) => assignment_operator

  procedure :: final => atlas_StructuredGrid_ComputeStencil__final
END TYPE

! Better not use, use setup member function instead !
!interface atlas_StructuredGrid_ComputeStencil
!  module procedure atlas_StructuredGrid_ComputeStencil__ctor
!end interface


!========================================================
contains
!========================================================


! -----------------------------------------------------------------------------
! Destructor

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_StructuredGrid_ComputeNorth__final_auto(this)
  type(atlas_StructuredGrid_ComputeNorth), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_StructuredGrid_ComputeWest__final_auto(this)
  type(atlas_StructuredGrid_ComputeWest), intent(inout) :: this
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif


! -----------------------------------------------------------------------------
! Constructors

function atlas_StructuredGrid_ComputeNorth__ctor(grid, halo) result(this)
  use, intrinsic :: iso_c_binding, only : c_int
  use fckit_c_interop_module, only: c_str
  use atlas_grid_StencilComputer_c_binding, only : atlas__grid__ComputeNorth__new, atlas__grid__ComputeNorth__delete
  use atlas_grid_module, only : atlas_StructuredGrid
  implicit none
  type(atlas_StructuredGrid_ComputeNorth) :: this
  type(atlas_StructuredGrid), intent(in) :: grid
  integer(c_int), intent(in) :: halo
  call this%reset_c_ptr( atlas__grid__ComputeNorth__new(grid%CPTR_PGIBUG_B, halo), &
    & fckit_c_deleter(atlas__grid__ComputeNorth__delete) )
  call this%return()
end function

function atlas_StructuredGrid_ComputeWest__ctor(grid, halo) result(this)
  use, intrinsic :: iso_c_binding, only : c_int
  use fckit_c_interop_module, only: c_str
  use atlas_grid_StencilComputer_c_binding, only : atlas__grid__ComputeWest__new, atlas__grid__ComputeWest__delete
  use atlas_grid_module, only : atlas_StructuredGrid
  implicit none
  type(atlas_StructuredGrid_ComputeWest) :: this
  type(atlas_StructuredGrid), intent(in) :: grid
  integer(c_int), intent(in) :: halo
  call this%reset_c_ptr( atlas__grid__ComputeWest__new(grid%CPTR_PGIBUG_B, halo), &
    & fckit_c_deleter(atlas__grid__ComputeWest__delete) )
  call this%return()
end function

subroutine atlas_StructuredGrid_ComputeStencil__setup(this, grid, stencil_width)
  use, intrinsic :: iso_c_binding, only : c_double
  use atlas_grid_module, only : atlas_StructuredGrid
  implicit none
  class(atlas_StructuredGrid_ComputeStencil) :: this
  type(atlas_StructuredGrid), intent(in) :: grid
  integer(ATLAS_KIND_IDX), intent(in) :: stencil_width
  this%stencil_width = stencil_width
  this%halo = (stencil_width + 1) / 2
  this%stencil_offset = stencil_width - floor(real(stencil_width,c_double) / 2._c_double + 1._c_double, ATLAS_KIND_IDX)
  this%compute_north = atlas_StructuredGrid_ComputeNorth(grid, this%halo)
  this%compute_west  = atlas_StructuredGrid_ComputeWest(grid, this%halo)
end subroutine



! ----------------------------------------------------------------------------------------

function atlas_StructuredGrid_ComputeNorth__execute_real32(this, y) result(index)
  use, intrinsic :: iso_c_binding, only : c_float
  use atlas_grid_StencilComputer_c_binding, only : atlas__grid__ComputeNorth__execute_real32
  implicit none
  integer(ATLAS_KIND_IDX) :: index
  class(atlas_StructuredGrid_ComputeNorth), intent(in) :: this
  real(c_float), intent(in) :: y
  index = atlas__grid__ComputeNorth__execute_real32(this%CPTR_PGIBUG_B, y) + 1
end function

function atlas_StructuredGrid_ComputeNorth__execute_real64(this, y) result(index)
  use, intrinsic :: iso_c_binding, only : c_double
  use atlas_grid_StencilComputer_c_binding, only : atlas__grid__ComputeNorth__execute_real64
  implicit none
  integer(ATLAS_KIND_IDX) :: index
  class(atlas_StructuredGrid_ComputeNorth), intent(in) :: this
  real(c_double), intent(in) :: y
  index = atlas__grid__ComputeNorth__execute_real64(this%CPTR_PGIBUG_B, y) + 1
end function
! ----------------------------------------------------------------------------------------

function atlas_StructuredGrid_ComputeWest__execute_real32(this, x, j) result(index)
  use, intrinsic :: iso_c_binding, only : c_float
  use atlas_grid_StencilComputer_c_binding, only : atlas__grid__ComputeWest__execute_real32
  implicit none
  integer(ATLAS_KIND_IDX) :: index
  class(atlas_StructuredGrid_ComputeWest), intent(in) :: this
  real(c_float), intent(in) :: x
  integer(ATLAS_KIND_IDX), intent(in) :: j
  index = atlas__grid__ComputeWest__execute_real32(this%CPTR_PGIBUG_B, x, j-int(1,ATLAS_KIND_IDX)) + 1
end function

function atlas_StructuredGrid_ComputeWest__execute_real64(this, x, j) result(index)
  use, intrinsic :: iso_c_binding, only : c_double
  use atlas_grid_StencilComputer_c_binding, only : atlas__grid__ComputeWest__execute_real64
  implicit none
  integer(ATLAS_KIND_IDX) :: index
  class(atlas_StructuredGrid_ComputeWest), intent(in) :: this
  real(c_double), intent(in) :: x
  integer(ATLAS_KIND_IDX), intent(in) :: j
  index = atlas__grid__ComputeWest__execute_real64(this%CPTR_PGIBUG_B, x, j-int(1,ATLAS_KIND_IDX)) + 1
end function
! ----------------------------------------------------------------------------------------


subroutine atlas_StructuredGrid_ComputeStencil__execute_real64(this, x, y, stencil)
  use, intrinsic :: iso_c_binding, only : c_double
  implicit none
  class(atlas_StructuredGrid_ComputeStencil), intent(in) :: this
  real(c_double), intent(in) :: x, y
  type(atlas_StructuredGrid_Stencil), intent(inout) :: stencil

  integer(ATLAS_KIND_IDX) :: jj
  stencil%width = this%stencil_width

  stencil%j_begin = this%compute_north%execute(y) - this%stencil_offset
  do jj = 1_ATLAS_KIND_IDX, this%stencil_width
    stencil%i_begin(jj) = this%compute_west%execute(x, stencil%j_begin + jj - 1) - this%stencil_offset
  enddo
end subroutine
! ----------------------------------------------------------------------------------------

subroutine atlas_StructuredGrid_Stencil__write (stencil, unit, iotype, v_list, iostat, iomsg)
  implicit none
  class(atlas_StructuredGrid_Stencil), intent(in) :: stencil
    INTEGER, INTENT(IN) :: unit
    CHARACTER(*), INTENT(IN) :: iotype
    INTEGER, INTENT(IN)  :: v_list(:)
    INTEGER, INTENT(OUT) :: iostat
    CHARACTER(*), INTENT(INOUT) :: iomsg
    integer(ATLAS_KIND_IDX) :: jlat, jlon
    do jlat = 1, stencil%width
      write(unit,'(a,I0,a)',advance='no',IOSTAT=iostat) " j: ", stencil%j(jlat), "    i = "
      do jlon = 1, stencil%width
        write(unit,'(I3)',advance='no',IOSTAT=iostat) stencil%i(jlon,jlat)
      enddo
      write(0,'(a)',IOSTAT=iostat) new_line('a')
    enddo 
end subroutine

function atlas_StructuredGrid_Stencil__j(this, j_index) result(j)
  integer(ATLAS_KIND_IDX) :: j
  class(atlas_StructuredGrid_Stencil), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: j_index
  j = this%j_begin + (j_index-1)
end function

function atlas_StructuredGrid_Stencil__i(this, i_index, j_index) result(i)
  integer(ATLAS_KIND_IDX) :: i
  class(atlas_StructuredGrid_Stencil), intent(in) :: this
  integer(ATLAS_KIND_IDX) :: i_index
  integer(ATLAS_KIND_IDX) :: j_index
  i = this%i_begin(j_index) + (i_index-1)
end function

! ----------------------------------------------------------------------------------------

function atlas_StructuredGrid_ComputeStencil__ctor(grid, stencil_width) result(this)
  use atlas_grid_module, only : atlas_StructuredGrid
  implicit none
  type(atlas_StructuredGrid_ComputeStencil) :: this
  type(atlas_StructuredGrid), intent(in) :: grid
  integer(ATLAS_KIND_IDX), intent(in) :: stencil_width
  call this%setup(grid, stencil_width)
end function


subroutine atlas_StructuredGrid_ComputeStencil__assignment(this, other)
implicit none
class(atlas_StructuredGrid_ComputeStencil), intent(inout) :: this
class(atlas_StructuredGrid_ComputeStencil), intent(in) :: other
call this%final()
write(0,*) "owners = ", other%compute_north%owners()
this%compute_north = other%compute_north
this%compute_west = other%compute_west
this%stencil_width = other%stencil_width
this%stencil_offset = other%stencil_offset
this%halo = other%halo
end subroutine

subroutine atlas_StructuredGrid_ComputeStencil__final(this)
  class(atlas_StructuredGrid_ComputeStencil), intent(inout) :: this
  call this%compute_north%final()
  call this%compute_west%final()
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_StencilComputer_module
