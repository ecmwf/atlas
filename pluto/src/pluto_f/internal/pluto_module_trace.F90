! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_trace

implicit none
private

public :: pluto_trace_t

type :: pluto_trace_t
contains
    procedure, nopass :: enable  => pluto_trace_enable
    procedure, nopass :: enabled => pluto_trace_enabled
end type

contains

subroutine pluto_trace_enable(enable)
    use iso_c_binding, only: c_int
    logical, optional :: enable
    logical :: do_enable
    interface
        subroutine c_pluto_trace_enable(enable) bind(c)
            use iso_c_binding, only: c_int
            integer(c_int), value :: enable
        end subroutine
    end interface
    do_enable = .true.
    if (present(enable)) then
        do_enable = enable
    endif
    if (do_enable) then
        call c_pluto_trace_enable(1_c_int)
    else
        call c_pluto_trace_enable(0_c_int)
    endif
end subroutine

function pluto_trace_enabled()
    use iso_c_binding, only: c_int
    logical :: pluto_trace_enabled
    integer(c_int) :: enabled
    interface
        function c_pluto_trace_enabled() result(enabled) bind(c)
            use iso_c_binding, only: c_int
            integer(c_int) :: enabled
        end function
    end interface
    enabled = c_pluto_trace_enabled()
    if (enabled == 1) then
        pluto_trace_enabled = .true.
    else
        pluto_trace_enabled = .false.
    endif
end function

end module
