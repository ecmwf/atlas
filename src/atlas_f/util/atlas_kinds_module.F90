! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_kinds_module

use, intrinsic :: iso_c_binding, only: &
  & c_int, c_long, c_double, c_float

implicit none
private

private :: c_long, c_int, c_double, c_float

public :: ATLAS_KIND_GIDX
public :: ATLAS_KIND_IDX

public :: ATLAS_KIND_REAL64
public :: ATLAS_KIND_REAL32
public :: ATLAS_KIND_INT64
public :: ATLAS_KIND_INT32


#if ATLAS_BITS_GLOBAL == 32
integer, parameter :: ATLAS_KIND_GIDX = c_int
#elif ATLAS_BITS_GLOBAL == 64
integer, parameter :: ATLAS_KIND_GIDX = c_long
#else
#error ATLAS_BITS_GLOBAL must be either 32 or 64
#endif

#if ATLAS_BITS_LOCAL == 32
integer, parameter :: ATLAS_KIND_IDX = c_int
#elif ATLAS_BITS_LOCAL == 64
integer, parameter :: ATLAS_KIND_IDX = c_long
#else
#error ATLAS_BITS_LOCAL must be either 32 or 64
#endif

integer, parameter :: ATLAS_KIND_REAL64 = c_double
integer, parameter :: ATLAS_KIND_REAL32 = c_float
integer, parameter :: ATLAS_KIND_INT64  = c_long
integer, parameter :: ATLAS_KIND_INT32  = c_int

ENUM, bind(c)
  enumerator :: openmode
  enumerator :: app = 1
  enumerator :: out = 16
end ENUM

! =============================================================================
CONTAINS
! =============================================================================


! -----------------------------------------------------------------------------

end module atlas_kinds_module
