! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_Field

! Purpose :
! -------
!   *Field* : Object containing field data and Metadata

! Methods :
! -------
!   name : The name or tag this field was created with
!   data : Return the field as a fortran array of specified shape
!   Metadata : Return object that can contain a variety of Metadata

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: name => Field__name
  procedure :: datatype => Field__datatype
  procedure :: metadata => Field__metadata
  !procedure :: function_space => Field__function_space
  procedure, private :: shape_array => Field__shape_array
  procedure, private :: shape_idx   => Field__shape_idx
  procedure :: size => Field__size
  procedure :: rank => Field__rank
  procedure :: bytes => Field__bytes
  procedure :: levels => Field__levels
  procedure :: kind => Field__kind
  procedure, private :: access_data1_int32 => Field__access_data1_int32
  procedure, private :: access_data2_int32 => Field__access_data2_int32
  procedure, private :: access_data3_int32 => Field__access_data3_int32
  procedure, private :: access_data1_int64 => Field__access_data1_int64
  procedure, private :: access_data2_int64 => Field__access_data2_int64
  procedure, private :: access_data3_int64 => Field__access_data3_int64
  procedure, private :: access_data1_real32 => Field__access_data1_real32
  procedure, private :: access_data2_real32 => Field__access_data2_real32
  procedure, private :: access_data3_real32 => Field__access_data3_real32
  procedure, private :: access_data4_real32 => Field__access_data4_real32
  procedure, private :: access_data1_real64 => Field__access_data1_real64
  procedure, private :: access_data2_real64 => Field__access_data2_real64
  procedure, private :: access_data3_real64 => Field__access_data3_real64
  procedure, private :: access_data4_real64 => Field__access_data4_real64
  procedure, private :: access_data2_real64_bounds => Field__access_data2_real64_bounds
  procedure, private :: access_data3_real64_bounds => Field__access_data3_real64_bounds
  procedure, private :: access_data4_real64_bounds => Field__access_data4_real64_bounds
  generic :: shape => shape_array, shape_idx
  generic :: access_data => &
    & access_data1_int32, access_data1_int64, access_data1_real32, access_data1_real64, &
    & access_data2_int32, access_data2_int64, access_data2_real32, access_data2_real64, access_data2_real64_bounds, &
    & access_data3_int32, access_data3_int64, access_data3_real32, access_data3_real64, access_data3_real64_bounds, &
    &                                         access_data4_real32, access_data4_real64, access_data4_real64_bounds
END TYPE atlas_Field

interface atlas_Field
  module procedure atlas_Field__create
  module procedure atlas_Field__create_name_kind_shape
  module procedure atlas_Field__create_kind_shape
end interface

!------------------------------------------------------------------------------
