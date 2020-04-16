! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_GridDistribution_module

use atlas_Grid_module, only : atlas_Grid
use fckit_owned_object_module, only: fckit_owned_object
use, intrinsic :: iso_c_binding, only : c_ptr

implicit none

public :: atlas_GridDistribution

private

!-----------------------------
! atlas_GridDistribution     !
!-----------------------------


!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_GridDistribution

! Purpose :
! -------
!   *GridDistribution* : Object passed from atlas to inspect grid distribution

! Methods :
! -------

! Author :
! ------
!   12-Mar-2014 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: partition_int32 => GridDistribution__partition_int32
  procedure, private :: partition_int64 => GridDistribution__partition_int64
  procedure :: nb_partitions => GridDistribution__nb_partitions
  generic :: partition => partition_int32, partition_int64
#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_GridDistribution__final_auto
#endif
END TYPE atlas_GridDistribution

!------------------------------------------------------------------------------

interface atlas_GridDistribution
  module procedure atlas_GridDistribution__cptr
  module procedure atlas_GridDistribution__ctor
  module procedure atlas_GridDistribution__ctor_gridconfig
end interface

private :: c_ptr
private :: fckit_owned_object

!========================================================
contains
!========================================================
! -----------------------------------------------------------------------------
! GridDistribution routines

function atlas_GridDistribution__cptr( cptr ) result(this)
  use atlas_distribution_c_binding
  type(atlas_GridDistribution) :: this
  type(c_ptr) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_GridDistribution__ctor( part, part0 ) result(this)
  use atlas_distribution_c_binding
  use atlas_kinds_module, only : ATLAS_KIND_IDX
  type(atlas_GridDistribution) :: this
  integer, intent(in) :: part(:)
  integer, intent(in), optional :: part0
  integer(ATLAS_KIND_IDX) :: npts
  integer :: opt_part0
  opt_part0 = 0
  if( present(part0) ) opt_part0 = part0
  npts = size(part)
  call this%reset_c_ptr( atlas__GridDistribution__new(npts, part, opt_part0) )
  call this%return()
end function

function atlas_GridDistribution__ctor_gridconfig( grid, config ) result(this)
  use atlas_distribution_c_binding
  use atlas_Config_module, only : atlas_Config
  use atlas_Grid_module, only : atlas_ReducedGaussianGrid
  type(atlas_GridDistribution) :: this
  class(atlas_Grid), intent (in) :: grid
  type(atlas_Config), intent(in) :: config
  call this%reset_c_ptr( atlas__GridDistribution__new_gridconfig(grid%CPTR_PGIBUG_A, config%CPTR_PGIBUG_B) )
  call this%return()
end function

function GridDistribution__partition_int32(this, i) result(partition)
  use, intrinsic :: iso_c_binding, only: c_long, c_int
  use atlas_distribution_c_binding
  integer(c_int) :: partition
  class(atlas_GridDistribution), intent(in) :: this
  integer(c_int), intent(in) :: i
  partition = atlas__GridDistribution__partition_int32(this%CPTR_PGIBUG_A, i-1)
end function

function GridDistribution__partition_int64(this, i) result(partition)
  use, intrinsic :: iso_c_binding, only: c_long, c_int
  use atlas_distribution_c_binding
  integer(c_long) :: partition
  class(atlas_GridDistribution), intent(in) :: this
  integer(c_long), intent(in) :: i
  partition = atlas__GridDistribution__partition_int64(this%CPTR_PGIBUG_A, i-1)
end function

function GridDistribution__nb_partitions(this) result(nb_partitions)
  use, intrinsic :: iso_c_binding, only: c_long
  use atlas_distribution_c_binding
  class(atlas_GridDistribution), intent(in) :: this
  integer(c_long) :: nb_partitions
  nb_partitions = atlas__GridDistribution__nb_partitions (this%CPTR_PGIBUG_A)
end function

! ----------------------------------------------------------------------------------------

ATLAS_FINAL subroutine atlas_GridDistribution__final_auto(this)
  type(atlas_GridDistribution), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_GridDistribution__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine

! ----------------------------------------------------------------------------------------

end module atlas_GridDistribution_module
