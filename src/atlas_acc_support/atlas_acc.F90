module atlas_acc
use openacc
implicit none
private

public :: atlas_acc_get_num_devices
public :: atlas_acc_map_data
public :: atlas_acc_unmap_data
public :: atlas_acc_is_present
public :: atlas_acc_get_device_type
public :: atlas_acc_deviceptr

contains

function atlas_acc_compiler_id() bind(C,name="atlas_acc_compiler_id") result(compiler_id)
  use, intrinsic :: iso_c_binding, only : c_int
  integer(c_int) :: compiler_id
  ! compiler_id must match number in atlas_acc.h enum type
#ifdef _CRAYFTN
  compiler_id = 2 ! cray
#else
  compiler_id = 0 ! unknown
#endif
end function

function atlas_acc_get_num_devices() bind(C,name="atlas_acc_get_num_devices") result(num_devices)
  use, intrinsic :: iso_c_binding, only : c_int
  integer(c_int) :: num_devices
  integer(acc_device_kind) :: devicetype

  devicetype = acc_get_device_type()
  num_devices = acc_get_num_devices(devicetype)
end function

subroutine atlas_acc_map_data(data_arg, data_dev, bytes) bind(C,name="atlas_acc_map_data")
  use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
  type(*), dimension(*) :: data_arg
  type(c_ptr), value :: data_dev
  integer(c_size_t), value :: bytes
  call acc_map_data(data_arg, data_dev, bytes)
end subroutine

subroutine atlas_acc_unmap_data(data_arg) bind(C,name="atlas_acc_unmap_data")
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(*), dimension(*) :: data_arg
  call acc_unmap_data(data_arg)
end subroutine

function atlas_acc_is_present(data_arg, bytes) bind(C,name="atlas_acc_is_present") result(is_present)
  use, intrinsic :: iso_c_binding, only : c_size_t, c_ptr, c_char, c_int
  integer(c_int) :: is_present
  logical :: lpresent
  type(c_ptr), value :: data_arg
  integer(c_size_t), value :: bytes
  character(kind=c_char), pointer :: data_f(:)
  call c_f_pointer(data_arg, data_f,[bytes])
  lpresent = acc_is_present(data_f)
  is_present = 0
  if (lpresent) is_present = 1
end function

function atlas_acc_deviceptr(data_arg) bind(C,name="atlas_acc_deviceptr") result(deviceptr)
  use, intrinsic :: iso_c_binding, only : c_ptr
  type(*), dimension(*) :: data_arg
  type(c_ptr):: deviceptr
  deviceptr = acc_deviceptr(data_arg)
end function

function atlas_acc_get_device_type() bind(C,name="atlas_acc_get_device_type") result(devicetype)
  use, intrinsic :: iso_c_binding, only : c_int
  integer(c_int) :: devicetype
  integer(acc_device_kind) :: acc_devicetype
  acc_devicetype = acc_get_device_type()
  if (acc_devicetype == acc_device_host .or. acc_devicetype == acc_device_none) then
    devicetype = 0
  else 
    devicetype = 1
  endif
end function

end module 
