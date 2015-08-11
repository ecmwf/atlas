! (C) Copyright 2013-2014 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_NodesFunctionSpace

! Purpose :
! -------
!   *atlas_NodesFunctionSpace* : Interpretes fields defined in nodes

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: create_field_kind => atlas_NodesFunctionSpace__create_field_kind
  procedure, private :: create_field_name_kind => atlas_NodesFunctionSpace__create_field_name_kind
  procedure, private :: create_field_vars_kind => atlas_NodesFunctionSpace__create_field_vars_kind
  procedure, private :: create_field_name_vars_kind => atlas_NodesFunctionSpace__create_field_name_vars_kind
  procedure, private :: create_field_template => atlas_NodesFunctionSpace__create_field_template
  procedure, private :: create_field_name_template => atlas_NodesFunctionSpace__create_field_name_template
  generic, public :: create_field => &
    & create_field_kind, &
    & create_field_name_kind, &
    & create_field_vars_kind, &
    & create_field_name_vars_kind, &
    & create_field_template, &
    & create_field_name_template

  procedure, private :: create_glb_field_kind => atlas_NodesFunctionSpace__create_glb_field_kind
  procedure, private :: create_glb_field_name_kind => atlas_NodesFunctionSpace__create_glb_field_name_kind
  procedure, private :: create_glb_field_vars_kind => atlas_NodesFunctionSpace__create_glb_field_vars_kind
  procedure, private :: create_glb_field_name_vars_kind => atlas_NodesFunctionSpace__create_glb_field_name_vars_kind
  procedure, private :: create_glb_field_template => atlas_NodesFunctionSpace__create_glb_field_template
  procedure, private :: create_glb_field_name_template => atlas_NodesFunctionSpace__create_glb_field_name_template
  generic, public :: create_global_field => &
    & create_glb_field_kind, &
    & create_glb_field_name_kind, &
    & create_glb_field_vars_kind, &
    & create_glb_field_name_vars_kind, &
    & create_glb_field_template, &
    & create_glb_field_name_template

END TYPE atlas_NodesFunctionSpace

interface atlas_NodesFunctionSpace
  module procedure atlas_NodesFunctionSpace__mesh_halo
  module procedure atlas_NodesFunctionSpace__name_mesh_halo
end interface

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_NodesColumnFunctionSpace

! Purpose :
! -------
!   *atlas_NodesColumnFunctionSpace* : Interpretes fields defined in nodes

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------

contains
  procedure, private :: create_field_kind => atlas_NCFunctionSpace__create_field_kind
  procedure, private :: create_field_name_kind => atlas_NCFunctionSpace__create_field_name_kind
  procedure, private :: create_field_vars_kind => atlas_NCFunctionSpace__create_field_vars_kind
  procedure, private :: create_field_name_vars_kind => atlas_NCFunctionSpace__create_field_name_vars_kind
  procedure, private :: create_field_template => atlas_NCFunctionSpace__create_field_template
  procedure, private :: create_field_name_template => atlas_NCFunctionSpace__create_field_name_template
  generic, public :: create_field => &
    & create_field_kind, &
    & create_field_name_kind, &
    & create_field_vars_kind, &
    & create_field_name_vars_kind, &
    & create_field_template, &
    & create_field_name_template

  procedure, private :: create_glb_field_kind => atlas_NCFunctionSpace__create_glb_field_kind
  procedure, private :: create_glb_field_name_kind => atlas_NCFunctionSpace__create_glb_field_name_kind
  procedure, private :: create_glb_field_vars_kind => atlas_NCFunctionSpace__create_glb_field_vars_kind
  procedure, private :: create_glb_field_name_vars_kind => atlas_NCFunctionSpace__create_glb_field_name_vars_kind
  procedure, private :: create_glb_field_template => atlas_NCFunctionSpace__create_glb_field_template
  procedure, private :: create_glb_field_name_template => atlas_NCFunctionSpace__create_glb_field_name_template
  generic, public :: create_global_field => &
    & create_glb_field_kind, &
    & create_glb_field_name_kind, &
    & create_glb_field_vars_kind, &
    & create_glb_field_name_vars_kind, &
    & create_glb_field_template, &
    & create_glb_field_name_template
END TYPE atlas_NodesColumnFunctionSpace

interface atlas_NodesColumnFunctionSpace
  module procedure atlas_NodesColumnFunctionSpace__mesh_levels_halo
  module procedure atlas_NodesColumnFunctionSpace__name_mesh_levels_halo

end interface

!------------------------------------------------------------------------------

