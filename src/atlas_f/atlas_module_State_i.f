! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_State

! Purpose :
! -------
!   *atlas_State* :
!       Container type of fields that are defined on the same points
!       Describes how nodes are ordered
!       Describes how parallelisation for fields is done
!       Describes interpolation between nodes

! Methods :
! -------
!   name : The name or tag this function space was created with
!   create_field : Create a new real field in this function space with given name
!   remove : Remove a field with given name
!   field : Access to a field with given name
!   parallelise : Setup halo-exchange information
!   halo_exchange : Perform halo exchange on field_data

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

!-- Field
  procedure, public :: add    => atlas_State__add
  procedure, public :: remove => atlas_State__remove
  procedure, public :: has    => atlas_State__has
  procedure, public :: size    => atlas_State__size
  procedure, private :: field_by_name  => atlas_State__field_by_name
  procedure, private :: field_by_index => atlas_State__field_by_index
  generic, public :: field => field_by_name, field_by_index
  procedure, public :: metadata => atlas_State__metadata

END TYPE atlas_State

interface atlas_State
  module procedure atlas_State__new
  module procedure atlas_State__generate
end interface

!------------------------------------------------------------------------------
