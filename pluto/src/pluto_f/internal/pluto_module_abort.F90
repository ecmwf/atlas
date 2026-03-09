! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_abort

implicit none
private

public :: pluto_abort

contains

subroutine pluto_abort(message)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: message
    write(error_unit, *) "PLUTO ABORT: ", trim(message)
    error stop 1
end subroutine

end module
