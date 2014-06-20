! (C) Copyright 2013-2014 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: FunctionSpace_type

! Purpose :
! -------
!   *FunctionSpace* : 
!       Container type of fields that are defined on the same points
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done
!       Describes interpolation between nodes

! Methods :
! -------
!   name : The name or tag this function space was created with
!   create_field : Create a new real field in this function space with given name
!   remove_field : Remove a field with given name
!   field : Access to a field with given name
!   parallelise : Setup halo-exchange information
!   halo_exchange : Perform halo exchange on field_data

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: name => FunctionSpace__name
  procedure :: dof => FunctionSpace__dof
  procedure :: glb_dof => FunctionSpace__glb_dof
  procedure :: create_field => FunctionSpace__create_field
  procedure :: remove_field => FunctionSpace__remove_field
  procedure :: bounds => FunctionSpace__bounds
  procedure :: field => FunctionSpace__field
  procedure, public :: has_field => FunctionSpace__has_field
  procedure :: parallelise => FunctionSpace__parallelise
  procedure, private :: FunctionSpace__halo_exchange_int32_r1
  procedure, private :: FunctionSpace__halo_exchange_int32_r2
  procedure, private :: FunctionSpace__halo_exchange_int32_r3
  procedure, private :: FunctionSpace__halo_exchange_real32_r1
  procedure, private :: FunctionSpace__halo_exchange_real32_r2
  procedure, private :: FunctionSpace__halo_exchange_real32_r3
  procedure, private :: FunctionSpace__halo_exchange_real64_r1
  procedure, private :: FunctionSpace__halo_exchange_real64_r2
  procedure, private :: FunctionSpace__halo_exchange_real64_r3
  procedure, private :: FunctionSpace__halo_exchange_real64_r4
  procedure :: get_halo_exchange => FunctionSpace__get_halo_exchange
  procedure :: get_gather => FunctionSpace__get_gather
  generic :: halo_exchange => &
      & FunctionSpace__halo_exchange_int32_r1, &
      & FunctionSpace__halo_exchange_int32_r2, &
      & FunctionSpace__halo_exchange_int32_r3, &
      & FunctionSpace__halo_exchange_real32_r1, &
      & FunctionSpace__halo_exchange_real32_r2, &
      & FunctionSpace__halo_exchange_real32_r3, &
      & FunctionSpace__halo_exchange_real64_r1, &
      & FunctionSpace__halo_exchange_real64_r2, &
      & FunctionSpace__halo_exchange_real64_r3, &
      & FunctionSpace__halo_exchange_real64_r4
      procedure, private :: FunctionSpace__gather_real32_r1
      procedure, private :: FunctionSpace__gather_real32_r2
      procedure, private :: FunctionSpace__gather_real32_r3
  procedure, private :: FunctionSpace__gather_real64_r1
  procedure, private :: FunctionSpace__gather_real64_r2
  procedure, private :: FunctionSpace__gather_real64_r3
  procedure, private :: FunctionSpace__gather_int32_r1
  procedure, private :: FunctionSpace__gather_int32_r2
  generic :: gather => &
      & FunctionSpace__gather_real32_r1, &
      & FunctionSpace__gather_real32_r2, &
      & FunctionSpace__gather_real32_r3, &
      & FunctionSpace__gather_real64_r1, &
      & FunctionSpace__gather_real64_r2, &
      & FunctionSpace__gather_real64_r3, &
      & FunctionSpace__gather_int32_r1, &
      & FunctionSpace__gather_int32_r2
END TYPE FunctionSpace_type

interface new_FunctionSpace
  module procedure new_FunctionSpace
end interface
!------------------------------------------------------------------------------
