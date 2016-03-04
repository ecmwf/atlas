! (C) Copyright 2013-2015 ECMWF.

#include "atlas/atlas_f.h"

module sb_mod

type, public :: atlas_object
  integer,public :: cpp_object_ptr
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
res = v1%cpp_object_ptr + v2%cpp_object_ptr
end program sb_program

