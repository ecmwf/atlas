! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_object) :: atlas_NodesFunctionSpace

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

  procedure, private :: create_field_name_kind => atlas_NodesFunctionSpace__create_field_name_kind
  procedure, private :: create_field_name_kind_lev => atlas_NodesFunctionSpace__create_field_name_kind_lev
  procedure, private :: create_field_name_kind_vars => atlas_NodesFunctionSpace__create_field_name_kind_vars
  procedure, private :: create_field_name_kind_lev_vars => atlas_NodesFunctionSpace__create_field_name_kind_lev_vars
  procedure, private :: create_field_name_template => atlas_NodesFunctionSpace__create_field_name_template
  generic, public :: create_field => &
    & create_field_name_kind, &
    & create_field_name_kind_lev, &
    & create_field_name_kind_vars, &
    & create_field_name_kind_lev_vars, &
    & create_field_name_template

  procedure, private :: create_glb_field_name_kind => atlas_NodesFunctionSpace__create_glb_field_name_kind
  procedure, private :: create_glb_field_name_kind_lev => atlas_NodesFunctionSpace__create_glb_field_name_kind_lev
  procedure, private :: create_glb_field_name_kind_vars => atlas_NodesFunctionSpace__create_glb_field_name_kind_vars
  procedure, private :: create_glb_field_name_kind_lev_vars => atlas_NodesFunctionSpace__create_glb_field_name_kind_lev_vars
  procedure, private :: create_glb_field_name_template => atlas_NodesFunctionSpace__create_glb_field_name_template
  generic, public :: create_global_field => &
    & create_glb_field_name_kind, &
    & create_glb_field_name_kind_lev, &
    & create_glb_field_name_kind_vars, &
    & create_glb_field_name_kind_lev_vars, &
    & create_glb_field_name_template

  procedure, private :: halo_exchange_fieldset => atlas_NodesFunctionSpace__halo_exchange_fieldset
  procedure, private :: halo_exchange_field => atlas_NodesFunctionSpace__halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_fieldset, halo_exchange_field
  procedure, public :: get_halo_exchange => atlas_NodesFunctionSpace__get_halo_exchange

  procedure, private :: gather_fieldset => atlas_NodesFunctionSpace__gather_fieldset
  procedure, private :: gather_field => atlas_NodesFunctionSpace__gather_field
  generic, public :: gather => gather_fieldset, gather_field
  procedure, public :: get_gather => atlas_NodesFunctionSpace__get_gather

  procedure, private :: scatter_fieldset => atlas_NodesFunctionSpace__scatter_fieldset
  procedure, private :: scatter_field => atlas_NodesFunctionSpace__scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field
  procedure, public :: get_scatter => atlas_NodesFunctionSpace__get_scatter

  procedure, private :: checksum_fieldset => atlas_NodesFunctionSpace__checksum_fieldset
  procedure, private :: checksum_field => atlas_NodesFunctionSpace__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field
  procedure, public :: get_checksum => atlas_NodesFunctionSpace__get_checksum

  procedure, private :: sum_real64_r0 => atlas_NodesFunctionSpace__sum_real64_r0
  procedure, private :: sum_real32_r0 => atlas_NodesFunctionSpace__sum_real32_r0
  procedure, private :: sum_int64_r0  => atlas_NodesFunctionSpace__sum_int64_r0
  procedure, private :: sum_int32_r0  => atlas_NodesFunctionSpace__sum_int32_r0
  procedure, private :: sum_real64_r1 => atlas_NodesFunctionSpace__sum_real64_r1
  procedure, private :: sum_real32_r1 => atlas_NodesFunctionSpace__sum_real32_r1
  procedure, private :: sum_int64_r1  => atlas_NodesFunctionSpace__sum_int64_r1
  procedure, private :: sum_int32_r1  => atlas_NodesFunctionSpace__sum_int32_r1
  procedure, private :: order_independent_sum_real32_r0 => atlas_NodesFunctionSpace__order_independent_sum_real32_r0
  procedure, private :: order_independent_sum_real64_r0 => atlas_NodesFunctionSpace__order_independent_sum_real64_r0
  procedure, private :: order_independent_sum_real32_r1 => atlas_NodesFunctionSpace__order_independent_sum_real32_r1
  procedure, private :: order_independent_sum_real64_r1 => atlas_NodesFunctionSpace__order_independent_sum_real64_r1
  procedure, private :: minimum_real32_r0 => atlas_NodesFunctionSpace__minimum_real32_r0
  procedure, private :: minimum_real64_r0 => atlas_NodesFunctionSpace__minimum_real64_r0
  procedure, private :: minimum_int32_r0  => atlas_NodesFunctionSpace__minimum_int32_r0
  procedure, private :: minimum_int64_r0  => atlas_NodesFunctionSpace__minimum_int64_r0
  procedure, private :: minimum_real32_r1 => atlas_NodesFunctionSpace__minimum_real32_r1
  procedure, private :: minimum_real64_r1 => atlas_NodesFunctionSpace__minimum_real64_r1
  procedure, private :: minimum_int32_r1  => atlas_NodesFunctionSpace__minimum_int32_r1
  procedure, private :: minimum_int64_r1  => atlas_NodesFunctionSpace__minimum_int64_r1
  procedure, private :: maximum_real32_r0 => atlas_NodesFunctionSpace__maximum_real32_r0
  procedure, private :: maximum_real64_r0 => atlas_NodesFunctionSpace__maximum_real64_r0
  procedure, private :: maximum_int32_r0  => atlas_NodesFunctionSpace__maximum_int32_r0
  procedure, private :: maximum_int64_r0  => atlas_NodesFunctionSpace__maximum_int64_r0
  procedure, private :: maximum_real32_r1 => atlas_NodesFunctionSpace__maximum_real32_r1
  procedure, private :: maximum_real64_r1 => atlas_NodesFunctionSpace__maximum_real64_r1
  procedure, private :: maximum_int32_r1  => atlas_NodesFunctionSpace__maximum_int32_r1
  procedure, private :: maximum_int64_r1  => atlas_NodesFunctionSpace__maximum_int64_r1
  procedure, private :: minimum_and_location_real32_r0 => atlas_NodesFunctionSpace__minloclev_real32_r0
  procedure, private :: minimum_and_location_real64_r0 => atlas_NodesFunctionSpace__minloclev_real64_r0
  procedure, private :: minimum_and_location_int32_r0  => atlas_NodesFunctionSpace__minloclev_int32_r0
  procedure, private :: minimum_and_location_int64_r0  => atlas_NodesFunctionSpace__minloclev_int64_r0
  procedure, private :: minimum_and_location_real32_r1 => atlas_NodesFunctionSpace__minloclev_real32_r1
  procedure, private :: minimum_and_location_real64_r1 => atlas_NodesFunctionSpace__minloclev_real64_r1
  procedure, private :: minimum_and_location_int32_r1  => atlas_NodesFunctionSpace__minloclev_int32_r1
  procedure, private :: minimum_and_location_int64_r1  => atlas_NodesFunctionSpace__minloclev_int64_r1
  procedure, private :: maximum_and_location_real32_r0 => atlas_NodesFunctionSpace__maxloclev_real32_r0
  procedure, private :: maximum_and_location_real64_r0 => atlas_NodesFunctionSpace__maxloclev_real64_r0
  procedure, private :: maximum_and_location_int32_r0  => atlas_NodesFunctionSpace__maxloclev_int32_r0
  procedure, private :: maximum_and_location_int64_r0  => atlas_NodesFunctionSpace__maxloclev_int64_r0
  procedure, private :: maximum_and_location_real32_r1 => atlas_NodesFunctionSpace__maxloclev_real32_r1
  procedure, private :: maximum_and_location_real64_r1 => atlas_NodesFunctionSpace__maxloclev_real64_r1
  procedure, private :: maximum_and_location_int32_r1  => atlas_NodesFunctionSpace__maxloclev_int32_r1
  procedure, private :: maximum_and_location_int64_r1  => atlas_NodesFunctionSpace__maxloclev_int64_r1
  procedure, private :: mean_real32_r0 => atlas_NodesFunctionSpace__mean_real32_r0
  procedure, private :: mean_real64_r0 => atlas_NodesFunctionSpace__mean_real64_r0
  procedure, private :: mean_int32_r0  => atlas_NodesFunctionSpace__mean_int32_r0
  procedure, private :: mean_int64_r0  => atlas_NodesFunctionSpace__mean_int64_r0
  procedure, private :: mean_real32_r1 => atlas_NodesFunctionSpace__mean_real32_r1
  procedure, private :: mean_real64_r1 => atlas_NodesFunctionSpace__mean_real64_r1
  procedure, private :: mean_int32_r1  => atlas_NodesFunctionSpace__mean_int32_r1
  procedure, private :: mean_int64_r1  => atlas_NodesFunctionSpace__mean_int64_r1
  procedure, private :: mean_and_stddev_real32_r0 => atlas_NodesFunctionSpace__mean_and_stddev_real32_r0
  procedure, private :: mean_and_stddev_real64_r0 => atlas_NodesFunctionSpace__mean_and_stddev_real64_r0
  procedure, private :: mean_and_stddev_int32_r0  => atlas_NodesFunctionSpace__mean_and_stddev_int32_r0
  procedure, private :: mean_and_stddev_int64_r0  => atlas_NodesFunctionSpace__mean_and_stddev_int64_r0
  procedure, private :: mean_and_stddev_real32_r1 => atlas_NodesFunctionSpace__mean_and_stddev_real32_r1
  procedure, private :: mean_and_stddev_real64_r1 => atlas_NodesFunctionSpace__mean_and_stddev_real64_r1
  procedure, private :: mean_and_stddev_int32_r1  => atlas_NodesFunctionSpace__mean_and_stddev_int32_r1
  procedure, private :: mean_and_stddev_int64_r1  => atlas_NodesFunctionSpace__mean_and_stddev_int64_r1

  generic, public :: minimum => &
    & minimum_real32_r0, minimum_real32_r1, &
    & minimum_real64_r0, minimum_real64_r1, &
    & minimum_int32_r0,  minimum_int32_r1,  &
    & minimum_int64_r0,  minimum_int64_r1

  procedure, public :: minimum_per_level => atlas_NodesFunctionSpace__minimum_per_level

  generic, public :: maximum => &
    & maximum_real32_r0, maximum_real32_r1, &
    & maximum_real64_r0, maximum_real64_r1, &
    & maximum_int32_r0,  maximum_int32_r1,  &
    & maximum_int64_r0,  maximum_int64_r1

  procedure, public :: maximum_per_level => atlas_NodesFunctionSpace__maximum_per_level

  generic, public :: minimum_and_location => &
    & minimum_and_location_real32_r0, minimum_and_location_real32_r1, &
    & minimum_and_location_real64_r0, minimum_and_location_real64_r1, &
    & minimum_and_location_int32_r0,  minimum_and_location_int32_r1,  &
    & minimum_and_location_int64_r0,  minimum_and_location_int64_r1

  procedure, public :: minimum_and_location_per_level => atlas_NodesFunctionSpace__minloc_per_level

  generic, public :: maximum_and_location => &
    & maximum_and_location_real32_r0, maximum_and_location_real32_r1, &
    & maximum_and_location_real64_r0, maximum_and_location_real64_r1, &
    & maximum_and_location_int32_r0,  maximum_and_location_int32_r1,  &
    & maximum_and_location_int64_r0,  maximum_and_location_int64_r1

  procedure, public :: maximum_and_location_per_level => atlas_NodesFunctionSpace__maxloc_per_level

  generic, public :: sum => &
    & sum_real32_r0, sum_real32_r1, &
    & sum_real64_r0, sum_real64_r1, &
    & sum_int32_r0,  sum_int32_r1,  &
    & sum_int64_r0,  sum_int64_r1

  procedure, public :: sum_per_level => atlas_NodesFunctionSpace__sum_per_level

  generic, public :: order_independent_sum => &
    & order_independent_sum_real32_r0, order_independent_sum_real32_r1, &
    & order_independent_sum_real64_r0, order_independent_sum_real64_r1


  procedure, public :: order_independent_sum_per_level => atlas_NodesFunctionSpace__order_independent_sum_per_level

  generic, public :: mean => &
    & mean_real32_r0, mean_real32_r1, &
    & mean_real64_r0, mean_real64_r1, &
    & mean_int32_r0,  mean_int32_r1,  &
    & mean_int64_r0,  mean_int64_r1

  procedure, public :: mean_per_level => atlas_NodesFunctionSpace__mean_per_level

  generic, public :: mean_and_standard_deviation => &
    & mean_and_stddev_real32_r0, mean_and_stddev_real32_r1, &
    & mean_and_stddev_real64_r0, mean_and_stddev_real64_r1, &
    & mean_and_stddev_int32_r0,  mean_and_stddev_int32_r1,  &
    & mean_and_stddev_int64_r0,  mean_and_stddev_int64_r1

  procedure, public :: mean_and_standard_deviation_per_level => atlas_NodesFunctionSpace__mean_and_stddev_per_level

  procedure :: finalize => atlas_NodesFunctionSpace__finalize
#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_NodesFunctionSpace__final
#endif


END TYPE atlas_NodesFunctionSpace

interface atlas_NodesFunctionSpace
  module procedure atlas_NodesFunctionSpace__cptr
  module procedure atlas_NodesFunctionSpace__mesh_halo
  module procedure atlas_NodesFunctionSpace__name_mesh_halo
end interface

interface assignment(=)
  module procedure atlas_NodesFunctionSpace__reset
end interface


!------------------------------------------------------------------------------

