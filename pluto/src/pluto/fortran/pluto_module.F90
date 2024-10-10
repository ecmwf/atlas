module pluto_module

use iso_c_binding, only: c_loc, c_ptr, c_int
implicit none
private

public :: pluto

interface
    subroutine c_pluto_host_set_default_resource(chars, n) bind(c)
        !use iso_c_binding, only: c_ptr, c_int
        type(c_ptr), value, intent(in) :: chars
        integer(c_int), intent(in) :: n
    end subroutine
end interface

type pluto_host_t
contains
    !procedure, nopass :: set_default_resource -> pluto_host_set_default_resource
end type

type pluto_t
    type(pluto_host_t) :: host
end type

contains

    subroutine pluto_host_set_default_resource(chars)
        character(len=*), target, intent(in) :: chars
        call c_pluto_host_set_default_resource(c_loc(chars), len(chars,kind=c_int))
    end subroutine

end module
