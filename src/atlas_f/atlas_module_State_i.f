! (C) Copyright 2013-2014 ECMWF.

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
!   remove_field : Remove a field with given name
!   field : Access to a field with given name
!   parallelise : Setup halo-exchange information
!   halo_exchange : Perform halo exchange on field_data

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

!-- Field
  procedure, public :: add    => atlas_State__add_field
  procedure, public :: remove_field => atlas_State__remove_field
  procedure, public :: has_field    => atlas_State__has_field
  procedure, public :: nb_fields    => atlas_State__nb_fields
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
