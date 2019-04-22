! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_JSON_module

implicit none

public :: atlas_PathName
public :: atlas_JSON

private


TYPE :: atlas_PathName
  character(len=1), allocatable, private :: string(:)
contains
  procedure :: str => atlas_PathName__str
END TYPE atlas_PathName

interface atlas_PathName
  module procedure atlas_PathName__ctor_str
end interface



TYPE :: atlas_JSON
  character(len=1), allocatable, private :: string(:)
contains
  procedure :: str => atlas_JSON__str

  procedure, public :: delete => atlas_JSON__delete

END TYPE atlas_JSON

interface atlas_JSON
  module procedure atlas_JSON__ctor_str
  module procedure atlas_JSON__ctor_path
end interface

!------------------------------------------------------------------------------
!========================================================
contains
!========================================================


function atlas_PathName__ctor_str(str) result(PathName)
  type(atlas_PathName) :: PathName
  character(len=*), intent(in) :: str
  integer i, nchars
  nchars = len(str)
  allocate( PathName%string(nchars) )
  do i=1,nchars
    PathName%string(i) = str(i:i)
  enddo
end function

function atlas_PathName__str(this) result(str)
  character(len=:), allocatable :: str
  class(atlas_PathName) :: this
  integer i, nchars
  nchars = size(this%string)
  allocate(character(len=nchars) :: str)
  do i=1,nchars
    str(i:i) = this%string(i)
  enddo
end function


function atlas_JSON__ctor_str(str) result(JSON)
  type(atlas_JSON) :: JSON
  character(len=*), intent(in) :: str
  integer i, nchars
  nchars = len(str)
  allocate( JSON%string(nchars) )
  do i=1,nchars
    JSON%string(i) = str(i:i)
  enddo
end function

function atlas_JSON__ctor_path(path) result(JSON)
  use, intrinsic :: iso_c_binding, only : c_char, c_int, c_ptr
  use atlas_atlas_read_file_c_binding
  use fckit_c_interop_module, only : c_str, c_ptr_to_string, c_ptr_free
  type(atlas_JSON) :: JSON
  type(atlas_PathName), intent(in) :: path
  character(len=:), allocatable :: str
  integer (c_int) :: iret
  type(c_ptr) :: str_cptr
  integer(c_int) :: str_size
  iret = atlas__read_file(c_str(path%str()), str_cptr, str_size)
  allocate(character(len=str_size) :: str )
  str = c_ptr_to_string(str_cptr)
  call c_ptr_free(str_cptr)
  JSON = atlas_JSON__ctor_str(str)
end function

function atlas_JSON__str(this) result(str)
  character(len=:), allocatable :: str
  class(atlas_JSON) :: this
  integer i, nchars
  nchars = size(this%string)
  allocate(character(len=nchars) :: str)
  do i=1,nchars
    str(i:i) = this%string(i)
  enddo
end function

subroutine atlas_JSON__delete(this)
  class(atlas_JSON), intent(inout) :: this
  if( allocated(this%string) ) deallocate(this%string)
end subroutine

end module atlas_JSON_module

