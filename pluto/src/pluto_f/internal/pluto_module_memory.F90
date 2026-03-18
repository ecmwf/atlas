! (C) Copyright 2016- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module pluto_module_memory

implicit none
private

public :: pluto_memory_t

type pluto_memory_t
contains
    procedure, nopass :: report => pluto_memory_report
end type

contains

function c_ptr_to_string(str_c_ptr,str_size) result(string)
  use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_size_t, c_f_pointer
  type(c_ptr), intent(in) :: str_c_ptr
  integer(c_size_t), intent(in) :: str_size
  character(kind=c_char,len=:), allocatable :: string
  character(kind=c_char,len=1), pointer  :: str_f_ptr(:)
  integer :: c
  call c_f_pointer( str_c_ptr , str_f_ptr, [str_size] )
  allocate( character(len=(str_size)) :: string )
  do c=1,str_size
    string(c:c) = str_f_ptr(c)
  enddo
end function

function pluto_memory_report() result(string)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    character(len=:), allocatable :: string
    interface
        subroutine c_pluto_memory_report(str_c_ptr, str_size) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), intent(out) :: str_c_ptr
            integer(c_size_t), intent(out) :: str_size
        end subroutine
        subroutine c_pluto_str_delete(str_c_ptr, str_size) bind(c)
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            type(c_ptr), value, intent(in) :: str_c_ptr
            integer(c_size_t), value, intent(in) :: str_size
        end subroutine
    end interface
    type(c_ptr) :: str_c_ptr
    integer(c_size_t) :: str_size
    call c_pluto_memory_report(str_c_ptr, str_size)
    string = c_ptr_to_string(str_c_ptr, str_size)
    call c_pluto_str_delete(str_c_ptr, str_size)
end function

end module
