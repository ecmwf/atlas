! (C) Copyright 2013-2014 ECMWF.

subroutine do_something_with_T2(v2)
    type(T2) :: v2
    v2%cpp_object_ptr = 2
end subroutine
