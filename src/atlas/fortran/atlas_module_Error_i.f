! (C) Copyright 2013-2014 ECMWF.

TYPE :: CodeLocation_type
  integer :: line
  character(len=:), allocatable :: file
  character(len=:), allocatable :: function
contains
ENDTYPE

interface atlas_abort
  module procedure atlas_abort_null
  module procedure atlas_abort_msg
  module procedure atlas_abort_msg_loc
end interface atlas_abort

interface atlas_throw_exception
  module procedure atlas_throw_exception_msg
  module procedure atlas_throw_exception_msg_loc
end interface atlas_throw_exception

interface atlas_throw_notimplemented
  module procedure atlas_throw_notimplemented_msg
  module procedure atlas_throw_notimplemented_msg_loc
  module procedure atlas_throw_notimplemented_loc
end interface atlas_throw_notimplemented

interface code_location
  module procedure code_location_null
  module procedure code_location_file_line
  module procedure code_location_file_line_func
end interface code_location

! Error codes
integer, parameter, public ::      &
  atlas_err_cleared         = -1,  &
  atlas_err_noerr           = 0,   &
  atlas_err_exception       = 1,   &
  atlas_err_usererror       = 2,   &
  atlas_err_seriousbug      = 3,   &
  atlas_err_notimplemented  = 4,   &
  atlas_err_assertionfailed = 5,   &
  atlas_err_badparameter    = 6,   &
  atlas_err_outofrange      = 7,   &
  atlas_err_stop            = 100, &
  atlas_err_abort           = 101, &
  atlas_err_cancel          = 102, &
  atlas_err_readerror       = 200, &
  atlas_err_writeerror      = 201
