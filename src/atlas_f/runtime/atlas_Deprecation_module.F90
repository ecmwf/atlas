! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module atlas_Deprecation_module

implicit none

public :: atlas_deprecation_errors
public :: atlas_deprecation_warnings

private

contains

logical function atlas_deprecation_warnings() result(enabled)
  logical, save :: initialized = .false.
  logical, save :: value = .true.
  character(len=16) :: env_value
  integer :: env_status
  integer :: env_length
  integer :: env_int
  integer :: read_status

  if (.not. initialized) then
    call get_environment_variable("ATLAS_DEPRECATION_WARNINGS", env_value, env_length, env_status)
    if (env_status == 0 .and. env_length > 0) then
      read(env_value(1:env_length), *, iostat=read_status) env_int
      if (read_status == 0) then
        value = (env_int /= 0)
      end if
    end if
    initialized = .true.
  end if

  enabled = value
end function

logical function atlas_deprecation_errors() result(enabled)
  logical, save :: initialized = .false.
  logical, save :: value = .false.
  character(len=16) :: env_value
  integer :: env_status
  integer :: env_length
  integer :: env_int
  integer :: read_status

  if (.not. initialized) then
    call get_environment_variable("ATLAS_DEPRECATION_ERRORS", env_value, env_length, env_status)
    if (env_status == 0 .and. env_length > 0) then
      read(env_value(1:env_length), *, iostat=read_status) env_int
      if (read_status == 0) then
        value = (env_int /= 0)
      end if
    end if
    initialized = .true.
  end if

  enabled = value
end function

end module atlas_Deprecation_module