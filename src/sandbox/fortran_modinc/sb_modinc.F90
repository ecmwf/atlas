! (C) Copyright 2013-2014 ECMWF.

#include "atlas/atlas_defines_fortran.h"


module sb_object_module

type, private :: private_type
  integer, public :: object
end type

type, public :: object_type
  type(private_type),public :: private
end type

end module


module sb_mod1
use sb_object_module
public T1
type, extends(object_type), public :: T1
end type
end module sb_mod1



module sb_mod2
use sb_object_module
public T2
type, extends(object_type), public :: T2
end type
end module sb_mod2


module sb_mod
use sb_mod1
use sb_mod2
end module sb_mod


program sb_program
use sb_mod
type(T1) :: v1
type(T2) :: v2
integer :: res
res = v1%private%object + v2%private%object
end program sb_program

