! (C) Copyright 2013-2015 ECMWF.

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

