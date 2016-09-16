subroutine atlas_write_to_fortran_unit(unit,msg_cptr) bind(C)
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr
  use fckit_c_interop_module, only : c_ptr_to_string
  integer(c_int), value, intent(in) :: unit
  type(c_ptr), value, intent(in) :: msg_cptr
  character(len=:), allocatable :: msg
  msg = c_ptr_to_string(msg_cptr)
  write(unit,'(A)',advance='no') msg
end subroutine
