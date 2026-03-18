! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_mpi

implicit none
private

public :: pluto_mpi_t

type pluto_mpi_t
contains
    procedure, nopass :: init => pluto_mpi_init
    procedure, nopass :: finalize => pluto_mpi_finalize
end type

contains

subroutine pluto_mpi_init()
    interface
        subroutine c_pluto_mpi_init() bind(c)
        end subroutine
    end interface
    call c_pluto_mpi_init()
end subroutine

subroutine pluto_mpi_finalize()
    interface
        subroutine c_pluto_mpi_finalize() bind(c)
        end subroutine
    end interface
    call c_pluto_mpi_finalize()
end subroutine

end module
