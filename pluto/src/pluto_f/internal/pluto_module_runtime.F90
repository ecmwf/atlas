! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_runtime

implicit none
private

public :: pluto_devices

contains

function pluto_devices()
    use, intrinsic :: iso_fortran_env, only: int32
    integer(int32) :: pluto_devices
    interface
        function c_pluto_devices() result(devices) bind(c)
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int) :: devices
        end function
    end interface
    pluto_devices = c_pluto_devices()
end function

end module
