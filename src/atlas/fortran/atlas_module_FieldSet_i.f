! (C) Copyright 2013-2014 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(object_type) :: FieldSet_type

! Purpose :
! -------
!   *FieldSet* : Object that groups Fields that go together
!       Fields can belong to several fieldsets simultaneously.
!       The actual ownership of the field lies in a FunctionSpace

! Methods :
! -------
!   add_field : The name or tag this field was created with
!   field : Return the field as a fortran array of specified shape
!   get_array : allocate a list of fields contained in the fieldset

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: size => FieldSet__size
  procedure, public :: add_field => FieldSet__add_field
  procedure, public :: has_field => FieldSet__has_field
  procedure, private :: field_by_name => FieldSet__field_by_name
  procedure, private :: field_by_idx => FieldSet__field_by_idx
  generic :: field => field_by_name, field_by_idx
  procedure, public :: get_array => FieldSet__fields
END TYPE FieldSet_type
!------------------------------------------------------------------------------
