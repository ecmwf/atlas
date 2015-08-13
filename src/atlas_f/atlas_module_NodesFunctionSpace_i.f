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

  procedure, public :: nb_nodes => atlas_NodesFunctionSpace__nb_nodes
  procedure, public :: mesh => atlas_NodesFunctionSpace__mesh
  procedure, public :: nodes => atlas_NodesFunctionSpace__nodes

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

  procedure, private :: halo_exchange_fieldset => atlas_NodesFunctionSpace__halo_exchange_fieldset
  procedure, private :: halo_exchange_field => atlas_NodesFunctionSpace__halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_fieldset, halo_exchange_field

  procedure, private :: gather_fieldset => atlas_NodesFunctionSpace__gather_fieldset
  procedure, private :: gather_field => atlas_NodesFunctionSpace__gather_field
  generic, public :: gather => gather_fieldset, gather_field

  procedure, private :: scatter_fieldset => atlas_NodesFunctionSpace__scatter_fieldset
  procedure, private :: scatter_field => atlas_NodesFunctionSpace__scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field

  procedure, private :: checksum_fieldset => atlas_NodesFunctionSpace__checksum_fieldset
  procedure, private :: checksum_field => atlas_NodesFunctionSpace__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field

  procedure, private :: sum_real64_r0 => atlas_NodesFunctionSpace__sum_real64_r0
  procedure, private :: sum_real32_r0 => atlas_NodesFunctionSpace__sum_real32_r0
  procedure, private :: sum_real64_r1 => atlas_NodesFunctionSpace__sum_real64_r1
  procedure, private :: sum_real32_r1 => atlas_NodesFunctionSpace__sum_real32_r1
  procedure, private :: order_independent_sum_real32_r0 => atlas_NodesFunctionSpace__order_independent_sum_real32_r0
  procedure, private :: order_independent_sum_real64_r0 => atlas_NodesFunctionSpace__order_independent_sum_real64_r0
  procedure, private :: order_independent_sum_real32_r1 => atlas_NodesFunctionSpace__order_independent_sum_real32_r1
  procedure, private :: order_independent_sum_real64_r1 => atlas_NodesFunctionSpace__order_independent_sum_real64_r1
  procedure, private :: minimum_real32_r0 => atlas_NodesFunctionSpace__minimum_real32_r0
  procedure, private :: minimum_real64_r0 => atlas_NodesFunctionSpace__minimum_real64_r0
  procedure, private :: minimum_real32_r1 => atlas_NodesFunctionSpace__minimum_real32_r1
  procedure, private :: minimum_real64_r1 => atlas_NodesFunctionSpace__minimum_real64_r1
  procedure, private :: maximum_real32_r0 => atlas_NodesFunctionSpace__maximum_real32_r0
  procedure, private :: maximum_real64_r0 => atlas_NodesFunctionSpace__maximum_real64_r0
  procedure, private :: maximum_real32_r1 => atlas_NodesFunctionSpace__maximum_real32_r1
  procedure, private :: maximum_real64_r1 => atlas_NodesFunctionSpace__maximum_real64_r1
  procedure, private :: minimum_and_location_real32_r0 => atlas_NodesFunctionSpace__minloc_real32_r0
  procedure, private :: minimum_and_location_real64_r0 => atlas_NodesFunctionSpace__minloc_real64_r0
  procedure, private :: minimum_and_location_real32_r1 => atlas_NodesFunctionSpace__minloc_real32_r1
  procedure, private :: minimum_and_location_real64_r1 => atlas_NodesFunctionSpace__minloc_real64_r1
  procedure, private :: maximum_and_location_real32_r0 => atlas_NodesFunctionSpace__maxloc_real32_r0
  procedure, private :: maximum_and_location_real64_r0 => atlas_NodesFunctionSpace__maxloc_real64_r0
  procedure, private :: maximum_and_location_real64_r1 => atlas_NodesFunctionSpace__maxloc_real64_r1
  procedure, private :: maximum_and_location_real32_r1 => atlas_NodesFunctionSpace__maxloc_real32_r1
  procedure, private :: mean_real32_r0 => atlas_NodesFunctionSpace__mean_real32_r0
  procedure, private :: mean_real64_r0 => atlas_NodesFunctionSpace__mean_real64_r0
  procedure, private :: mean_real32_r1 => atlas_NodesFunctionSpace__mean_real32_r1
  procedure, private :: mean_real64_r1 => atlas_NodesFunctionSpace__mean_real64_r1
  procedure, private :: mean_and_stddev_real32_r0 => atlas_NodesFunctionSpace__mean_and_stddev_real32_r0
  procedure, private :: mean_and_stddev_real64_r0 => atlas_NodesFunctionSpace__mean_and_stddev_real64_r0
  procedure, private :: mean_and_stddev_real32_r1 => atlas_NodesFunctionSpace__mean_and_stddev_real32_r1
  procedure, private :: mean_and_stddev_real64_r1 => atlas_NodesFunctionSpace__mean_and_stddev_real64_r1

  generic, public :: minimum => &
    & minimum_real32_r0, minimum_real32_r1, &
    & minimum_real64_r0, minimum_real64_r1

  generic, public :: maximum => &
    & maximum_real32_r0, maximum_real32_r1, &
    & maximum_real64_r0, maximum_real64_r1

  generic, public :: minimum_and_location => &
    & minimum_and_location_real32_r0, minimum_and_location_real32_r1, &
    & minimum_and_location_real64_r0, minimum_and_location_real64_r1

  generic, public :: maximum_and_location => &
    & maximum_and_location_real32_r0, maximum_and_location_real32_r1, &
    & maximum_and_location_real64_r0, maximum_and_location_real64_r1

  generic, public :: sum => &
    & sum_real32_r0, sum_real32_r1, &
    & sum_real64_r0, sum_real64_r1

  generic, public :: order_independent_sum => &
    & order_independent_sum_real32_r0, order_independent_sum_real32_r1, &
    & order_independent_sum_real64_r0, order_independent_sum_real64_r1

  generic, public :: mean => &
    & mean_real32_r0, mean_real32_r1, &
    & mean_real64_r0, mean_real64_r1

  generic, public :: mean_and_standard_deviation => &
    & mean_and_stddev_real32_r0, mean_and_stddev_real32_r1, &
    & mean_and_stddev_real64_r0, mean_and_stddev_real64_r1

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

  procedure, public :: nb_nodes => atlas_NodesColumnFunctionSpace__nb_nodes
  procedure, public :: nb_levels => atlas_NodesColumnFunctionSpace__nb_levels
  procedure, public :: mesh => atlas_NodesColumnFunctionSpace__mesh
  procedure, public :: nodes => atlas_NodesColumnFunctionSpace__nodes

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

  procedure, private :: halo_exchange_fieldset => atlas_NCFunctionSpace__halo_exchange_fieldset
  procedure, private :: halo_exchange_field => atlas_NCFunctionSpace__halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_fieldset, halo_exchange_field

  procedure, private :: gather_fieldset => atlas_NCFunctionSpace__gather_fieldset
  procedure, private :: gather_field => atlas_NCFunctionSpace__gather_field
  generic, public :: gather => gather_fieldset, gather_field

  procedure, private :: scatter_fieldset => atlas_NCFunctionSpace__scatter_fieldset
  procedure, private :: scatter_field => atlas_NCFunctionSpace__scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field

  procedure, private :: checksum_fieldset => atlas_NCFunctionSpace__checksum_fieldset
  procedure, private :: checksum_field => atlas_NCFunctionSpace__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field

END TYPE atlas_NodesColumnFunctionSpace

interface atlas_NodesColumnFunctionSpace
  module procedure atlas_NodesColumnFunctionSpace__mesh_levels_halo
  module procedure atlas_NodesColumnFunctionSpace__name_mesh_levels_halo

end interface

!------------------------------------------------------------------------------

