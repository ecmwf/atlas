! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------

TYPE, extends(object_type) :: atlas_PathName
  character(kind=c_char,len=1), allocatable, private :: string(:)
contains
  procedure :: str => atlas_PathName__str
END TYPE atlas_PathName

interface atlas_PathName
  module procedure atlas_PathName__ctor_str
end interface

!------------------------------------------------------------------------------
