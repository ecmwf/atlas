! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_scope

implicit none
private

public :: pluto_scope_t

type pluto_scope_t
contains
    procedure, nopass :: push => pluto_scope_push
    procedure, nopass :: pop  => pluto_scope_pop
end type

contains

subroutine pluto_scope_push()
    interface
        subroutine c_pluto_scope_push() bind(c)
        end subroutine
    end interface
    call c_pluto_scope_push()
end subroutine

subroutine pluto_scope_pop()
    interface
        subroutine c_pluto_scope_pop() bind(c)
        end subroutine
    end interface
    call c_pluto_scope_pop()
end subroutine

end module
