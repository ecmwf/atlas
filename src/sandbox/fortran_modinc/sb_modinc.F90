! (C) Copyright 2013-2014 ECMWF.

#include "atlas/atlas_defines_fortran.h"

module sb_mod

type, private :: private_type
  integer, public :: object
end type

type, public :: object_type
  type(private_type),public :: private
end type

#include "mod1.h"
#include "mod2.h"
contains
#include "mod1.f"
#include "mod2.f"
end module sb_mod


program sb_program
use sb_mod
type(T1) :: v1
type(T2) :: v2
integer :: res
res = v1%private%object + v2%private%object
end program sb_program

