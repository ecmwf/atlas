! (C) Copyright 2013-2014 ECMWF.

function atlas_PathName_str(str) result(PathName)
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

