! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_Metadata

! Purpose :
! -------
!   *Metadata* : Container of Metadata, parameters or attributes
!       The Metadata are added as key, value pairs

! Methods :
! -------
!   add : Add a new property with given key and value
!   set : Modify a property with given key and value
!   get : Return a property value for given key

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: add_logical => Metadata__add_logical
  procedure, private :: add_integer => Metadata__add_integer
  procedure, private :: add_real32 => Metadata__add_real32
  procedure, private :: add_real64 => Metadata__add_real64
  procedure, private :: add_string => Metadata__add_string
  procedure :: has => Metadata__has
  generic :: add => add_logical, add_integer, add_real32, add_real64, add_string
  generic :: set => add_logical, add_integer, add_real32, add_real64, add_string
  procedure :: get_integer => Metadata__get_integer
  procedure :: get_logical => Metadata__get_logical
  procedure :: get_real32 => Metadata__get_real32
  procedure :: get_real64 => Metadata__get_real64
  procedure :: get_string => Metadata__get_string
  generic :: get => get_integer, get_logical, get_real32, get_real64, get_string
END TYPE atlas_Metadata

!------------------------------------------------------------------------------

interface atlas_Metadata
  module procedure atlas_Metadata__ctor
end interface

!------------------------------------------------------------------------------
