! (C) Copyright 2013-2015 ECMWF.

#include "atlas/atlas_f.h"

module atlas_kinds_module

use, intrinsic :: iso_c_binding, only: &
  & c_int, c_long

implicit none
private

private :: c_long, c_int

public :: ATLAS_KIND_GIDX
public :: ATLAS_KIND_IDX



#if ATLAS_BITS_GLOBAL == 32
integer, parameter :: ATLAS_KIND_GIDX = c_int
#elif ATLAS_BITS_GLOBAL == 64
integer, parameter :: ATLAS_KIND_GIDX = c_long
#else
#error ATLAS_BITS_GLOBAL must be either 32 or 64
#endif

integer, parameter :: ATLAS_KIND_IDX = c_int

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
