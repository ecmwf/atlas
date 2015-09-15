! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_Value

! Purpose :
! -------
!   *Value* : Container of Value, which internally stores any standard type

! Methods :
! -------
!   get : Return a value in standard type

! Author :
! ------
!   June-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: get_logical => atlas_Value__get_logical
  procedure :: get_int32   => atlas_Value__get_int32
  procedure :: get_int64   => atlas_Value__get_int64
  procedure :: get_real32  => atlas_Value__get_real32
  procedure :: get_real64  => atlas_Value__get_real64
  procedure :: get_string  => atlas_Value__get_string
  procedure :: get_array_int32  => atlas_Value__get_array_int32
  procedure :: get_array_int64  => atlas_Value__get_array_int64
  procedure :: get_array_real32 => atlas_Value__get_array_real32
  procedure :: get_array_real64 => atlas_Value__get_array_real64
  generic :: get => get_logical, get_int32, get_int64, get_real32, get_real64, get_string, &
                    get_array_int32, get_array_int64, get_array_real32, get_array_real64

END TYPE atlas_Value

!------------------------------------------------------------------------------

interface atlas_Value
  module procedure atlas_Value__ctor_int32
  module procedure atlas_Value__ctor_int64
  module procedure atlas_Value__ctor_real32
  module procedure atlas_Value__ctor_real64
  module procedure atlas_Value__ctor_string
  module procedure atlas_Value__ctor_array_int32
  module procedure atlas_Value__ctor_array_int64
  module procedure atlas_Value__ctor_array_real32
  module procedure atlas_Value__ctor_array_real64
end interface

!------------------------------------------------------------------------------
