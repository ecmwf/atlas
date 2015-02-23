! (C) Copyright 2013-2014 ECMWF.

integer, parameter :: ATLAS_CODELOCATION_FILE_STRLEN     = 1024
integer, parameter :: ATLAS_CODELOCATION_FUNCTION_STRLEN = 1024

TYPE :: CodeLocation_type
  integer :: line
  character(len=ATLAS_CODELOCATION_FILE_STRLEN) :: file
  character(len=ATLAS_CODELOCATION_FILE_STRLEN) :: function
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
  module procedure atlas_throw_exception_loc
end interface atlas_throw_exception

interface atlas_throw_notimplemented
  module procedure atlas_throw_notimplemented_msg
  module procedure atlas_throw_notimplemented_loc
  module procedure atlas_throw_notimplemented_msg_loc
end interface atlas_throw_notimplemented

interface atlas_throw_outofrange
  module procedure atlas_throw_outofrange_msg
  module procedure atlas_throw_outofrange_loc
  module procedure atlas_throw_outofrange_msg_loc
  module procedure atlas_throw_outofrange_range
  module procedure atlas_throw_outofrange_range_loc
end interface atlas_throw_outofrange

interface atlas_throw_seriousbug
  module procedure atlas_throw_seriousbug_msg
  module procedure atlas_throw_seriousbug_loc
  module procedure atlas_throw_seriousbug_msg_loc
end interface atlas_throw_seriousbug

interface atlas_throw_usererror
  module procedure atlas_throw_usererror_msg
  module procedure atlas_throw_usererror_loc
  module procedure atlas_throw_usererror_msg_loc
end interface atlas_throw_usererror

interface atlas_throw_assertionfailed
  module procedure atlas_throw_assertionfailed_msg
  module procedure atlas_throw_assertionfailed_loc
  module procedure atlas_throw_assertionfailed_msg_loc
end interface atlas_throw_assertionfailed


interface code_location
  module procedure code_location_null
  module procedure code_location_file_line
  module procedure code_location_file_line_func
end interface code_location

! Error codes
integer, parameter, public ::      &
  atlas_err_cleared         =  1,  &
  atlas_err_noerr           =  0,   &
  atlas_err_exception       = -1,   &
  atlas_err_usererror       = -2,   &
  atlas_err_seriousbug      = -3,   &
  atlas_err_notimplemented  = -4,   &
  atlas_err_assertionfailed = -5,   &
  atlas_err_badparameter    = -6,   &
  atlas_err_outofrange      = -7,   &
  atlas_err_stop            = -100, &
  atlas_err_abort           = -101, &
  atlas_err_cancel          = -102, &
  atlas_err_readerror       = -200, &
  atlas_err_writeerror      = -201, &
  atlas_err_unknown         = -999
