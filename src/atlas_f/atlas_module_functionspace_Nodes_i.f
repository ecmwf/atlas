! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_NextFunctionSpace) :: atlas_functionspace_Nodes

! Purpose :
! -------
!   *atlas_functionspace_Nodes* : Interpretes fields defined in nodes

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains


  procedure, public :: nb_nodes => atlas_functionspace_Nodes__nb_nodes
  procedure, public :: mesh => atlas_functionspace_Nodes__mesh
  procedure, public :: nodes => atlas_functionspace_Nodes__nodes

  procedure, private :: create_field_name_kind => atlas_functionspace_Nodes__create_field_name_kind
  procedure, private :: create_field_name_kind_lev => atlas_functionspace_Nodes__create_field_name_kind_lev
  procedure, private :: create_field_name_kind_vars => atlas_functionspace_Nodes__create_field_name_kind_vars
  procedure, private :: create_field_name_kind_lev_vars => atlas_functionspace_Nodes__create_field_name_kind_lev_vars
  procedure, private :: create_field_name_template => atlas_functionspace_Nodes__create_field_name_template
  generic, public :: create_field => &
    & create_field_name_kind, &
    & create_field_name_kind_lev, &
    & create_field_name_kind_vars, &
    & create_field_name_kind_lev_vars, &
    & create_field_name_template

  procedure, private :: create_glb_field_name_kind => atlas_functionspace_Nodes__create_glb_field_name_kind
  procedure, private :: create_glb_field_name_kind_lev => atlas_functionspace_Nodes__create_glb_field_name_kind_lev
  procedure, private :: create_glb_field_name_kind_vars => atlas_functionspace_Nodes__create_glb_field_name_kind_vars
  procedure, private :: create_glb_field_name_kind_lev_vars => atlas_functionspace_Nodes__create_glb_field_name_kind_lev_vars
  procedure, private :: create_glb_field_name_template => atlas_functionspace_Nodes__create_glb_field_name_template
  generic, public :: create_global_field => &
    & create_glb_field_name_kind, &
    & create_glb_field_name_kind_lev, &
    & create_glb_field_name_kind_vars, &
    & create_glb_field_name_kind_lev_vars, &
    & create_glb_field_name_template

  procedure, private :: halo_exchange_fieldset => atlas_functionspace_Nodes__halo_exchange_fieldset
  procedure, private :: halo_exchange_field => atlas_functionspace_Nodes__halo_exchange_field
  generic, public :: halo_exchange => halo_exchange_fieldset, halo_exchange_field
  procedure, public :: get_halo_exchange => atlas_functionspace_Nodes__get_halo_exchange

  procedure, private :: gather_fieldset => atlas_functionspace_Nodes__gather_fieldset
  procedure, private :: gather_field => atlas_functionspace_Nodes__gather_field
  generic, public :: gather => gather_fieldset, gather_field
  procedure, public :: get_gather => atlas_functionspace_Nodes__get_gather

  procedure, private :: scatter_fieldset => atlas_functionspace_Nodes__scatter_fieldset
  procedure, private :: scatter_field => atlas_functionspace_Nodes__scatter_field
  generic, public :: scatter => scatter_fieldset, scatter_field
  procedure, public :: get_scatter => atlas_functionspace_Nodes__get_scatter

  procedure, private :: checksum_fieldset => atlas_functionspace_Nodes__checksum_fieldset
  procedure, private :: checksum_field => atlas_functionspace_Nodes__checksum_field
  generic, public :: checksum => checksum_fieldset, checksum_field
  procedure, public :: get_checksum => atlas_functionspace_Nodes__get_checksum

  procedure, private :: sum_real64_r0 => atlas_functionspace_Nodes__sum_real64_r0
  procedure, private :: sum_real32_r0 => atlas_functionspace_Nodes__sum_real32_r0
  procedure, private :: sum_int64_r0  => atlas_functionspace_Nodes__sum_int64_r0
  procedure, private :: sum_int32_r0  => atlas_functionspace_Nodes__sum_int32_r0
  procedure, private :: sum_real64_r1 => atlas_functionspace_Nodes__sum_real64_r1
  procedure, private :: sum_real32_r1 => atlas_functionspace_Nodes__sum_real32_r1
  procedure, private :: sum_int64_r1  => atlas_functionspace_Nodes__sum_int64_r1
  procedure, private :: sum_int32_r1  => atlas_functionspace_Nodes__sum_int32_r1
  procedure, private :: order_independent_sum_real32_r0 => atlas_functionspace_Nodes__order_independent_sum_real32_r0
  procedure, private :: order_independent_sum_real64_r0 => atlas_functionspace_Nodes__order_independent_sum_real64_r0
  procedure, private :: order_independent_sum_real32_r1 => atlas_functionspace_Nodes__order_independent_sum_real32_r1
  procedure, private :: order_independent_sum_real64_r1 => atlas_functionspace_Nodes__order_independent_sum_real64_r1
  procedure, private :: minimum_real32_r0 => atlas_functionspace_Nodes__minimum_real32_r0
  procedure, private :: minimum_real64_r0 => atlas_functionspace_Nodes__minimum_real64_r0
  procedure, private :: minimum_int32_r0  => atlas_functionspace_Nodes__minimum_int32_r0
  procedure, private :: minimum_int64_r0  => atlas_functionspace_Nodes__minimum_int64_r0
  procedure, private :: minimum_real32_r1 => atlas_functionspace_Nodes__minimum_real32_r1
  procedure, private :: minimum_real64_r1 => atlas_functionspace_Nodes__minimum_real64_r1
  procedure, private :: minimum_int32_r1  => atlas_functionspace_Nodes__minimum_int32_r1
  procedure, private :: minimum_int64_r1  => atlas_functionspace_Nodes__minimum_int64_r1
  procedure, private :: maximum_real32_r0 => atlas_functionspace_Nodes__maximum_real32_r0
  procedure, private :: maximum_real64_r0 => atlas_functionspace_Nodes__maximum_real64_r0
  procedure, private :: maximum_int32_r0  => atlas_functionspace_Nodes__maximum_int32_r0
  procedure, private :: maximum_int64_r0  => atlas_functionspace_Nodes__maximum_int64_r0
  procedure, private :: maximum_real32_r1 => atlas_functionspace_Nodes__maximum_real32_r1
  procedure, private :: maximum_real64_r1 => atlas_functionspace_Nodes__maximum_real64_r1
  procedure, private :: maximum_int32_r1  => atlas_functionspace_Nodes__maximum_int32_r1
  procedure, private :: maximum_int64_r1  => atlas_functionspace_Nodes__maximum_int64_r1
  procedure, private :: minimum_and_location_real32_r0 => atlas_functionspace_Nodes__minloclev_real32_r0
  procedure, private :: minimum_and_location_real64_r0 => atlas_functionspace_Nodes__minloclev_real64_r0
  procedure, private :: minimum_and_location_int32_r0  => atlas_functionspace_Nodes__minloclev_int32_r0
  procedure, private :: minimum_and_location_int64_r0  => atlas_functionspace_Nodes__minloclev_int64_r0
  procedure, private :: minimum_and_location_real32_r1 => atlas_functionspace_Nodes__minloclev_real32_r1
  procedure, private :: minimum_and_location_real64_r1 => atlas_functionspace_Nodes__minloclev_real64_r1
  procedure, private :: minimum_and_location_int32_r1  => atlas_functionspace_Nodes__minloclev_int32_r1
  procedure, private :: minimum_and_location_int64_r1  => atlas_functionspace_Nodes__minloclev_int64_r1
  procedure, private :: maximum_and_location_real32_r0 => atlas_functionspace_Nodes__maxloclev_real32_r0
  procedure, private :: maximum_and_location_real64_r0 => atlas_functionspace_Nodes__maxloclev_real64_r0
  procedure, private :: maximum_and_location_int32_r0  => atlas_functionspace_Nodes__maxloclev_int32_r0
  procedure, private :: maximum_and_location_int64_r0  => atlas_functionspace_Nodes__maxloclev_int64_r0
  procedure, private :: maximum_and_location_real32_r1 => atlas_functionspace_Nodes__maxloclev_real32_r1
  procedure, private :: maximum_and_location_real64_r1 => atlas_functionspace_Nodes__maxloclev_real64_r1
  procedure, private :: maximum_and_location_int32_r1  => atlas_functionspace_Nodes__maxloclev_int32_r1
  procedure, private :: maximum_and_location_int64_r1  => atlas_functionspace_Nodes__maxloclev_int64_r1
  procedure, private :: mean_real32_r0 => atlas_functionspace_Nodes__mean_real32_r0
  procedure, private :: mean_real64_r0 => atlas_functionspace_Nodes__mean_real64_r0
  procedure, private :: mean_int32_r0  => atlas_functionspace_Nodes__mean_int32_r0
  procedure, private :: mean_int64_r0  => atlas_functionspace_Nodes__mean_int64_r0
  procedure, private :: mean_real32_r1 => atlas_functionspace_Nodes__mean_real32_r1
  procedure, private :: mean_real64_r1 => atlas_functionspace_Nodes__mean_real64_r1
  procedure, private :: mean_int32_r1  => atlas_functionspace_Nodes__mean_int32_r1
  procedure, private :: mean_int64_r1  => atlas_functionspace_Nodes__mean_int64_r1
  procedure, private :: mean_and_stddev_real32_r0 => atlas_functionspace_Nodes__mean_and_stddev_real32_r0
  procedure, private :: mean_and_stddev_real64_r0 => atlas_functionspace_Nodes__mean_and_stddev_real64_r0
  procedure, private :: mean_and_stddev_int32_r0  => atlas_functionspace_Nodes__mean_and_stddev_int32_r0
  procedure, private :: mean_and_stddev_int64_r0  => atlas_functionspace_Nodes__mean_and_stddev_int64_r0
  procedure, private :: mean_and_stddev_real32_r1 => atlas_functionspace_Nodes__mean_and_stddev_real32_r1
  procedure, private :: mean_and_stddev_real64_r1 => atlas_functionspace_Nodes__mean_and_stddev_real64_r1
  procedure, private :: mean_and_stddev_int32_r1  => atlas_functionspace_Nodes__mean_and_stddev_int32_r1
  procedure, private :: mean_and_stddev_int64_r1  => atlas_functionspace_Nodes__mean_and_stddev_int64_r1

  generic, public :: minimum => &
    & minimum_real32_r0, minimum_real32_r1, &
    & minimum_real64_r0, minimum_real64_r1, &
    & minimum_int32_r0,  minimum_int32_r1,  &
    & minimum_int64_r0,  minimum_int64_r1

  procedure, public :: minimum_per_level => atlas_functionspace_Nodes__minimum_per_level

  generic, public :: maximum => &
    & maximum_real32_r0, maximum_real32_r1, &
    & maximum_real64_r0, maximum_real64_r1, &
    & maximum_int32_r0,  maximum_int32_r1,  &
    & maximum_int64_r0,  maximum_int64_r1

  procedure, public :: maximum_per_level => atlas_functionspace_Nodes__maximum_per_level

  generic, public :: minimum_and_location => &
    & minimum_and_location_real32_r0, minimum_and_location_real32_r1, &
    & minimum_and_location_real64_r0, minimum_and_location_real64_r1, &
    & minimum_and_location_int32_r0,  minimum_and_location_int32_r1,  &
    & minimum_and_location_int64_r0,  minimum_and_location_int64_r1

  procedure, public :: minimum_and_location_per_level => atlas_functionspace_Nodes__minloc_per_level

  generic, public :: maximum_and_location => &
    & maximum_and_location_real32_r0, maximum_and_location_real32_r1, &
    & maximum_and_location_real64_r0, maximum_and_location_real64_r1, &
    & maximum_and_location_int32_r0,  maximum_and_location_int32_r1,  &
    & maximum_and_location_int64_r0,  maximum_and_location_int64_r1

  procedure, public :: maximum_and_location_per_level => atlas_functionspace_Nodes__maxloc_per_level

  generic, public :: sum => &
    & sum_real32_r0, sum_real32_r1, &
    & sum_real64_r0, sum_real64_r1, &
    & sum_int32_r0,  sum_int32_r1,  &
    & sum_int64_r0,  sum_int64_r1

  procedure, public :: sum_per_level => atlas_functionspace_Nodes__sum_per_level

  generic, public :: order_independent_sum => &
    & order_independent_sum_real32_r0, order_independent_sum_real32_r1, &
    & order_independent_sum_real64_r0, order_independent_sum_real64_r1


  procedure, public :: order_independent_sum_per_level => atlas_functionspace_Nodes__order_independent_sum_per_level

  generic, public :: mean => &
    & mean_real32_r0, mean_real32_r1, &
    & mean_real64_r0, mean_real64_r1, &
    & mean_int32_r0,  mean_int32_r1,  &
    & mean_int64_r0,  mean_int64_r1

  procedure, public :: mean_per_level => atlas_functionspace_Nodes__mean_per_level

  generic, public :: mean_and_standard_deviation => &
    & mean_and_stddev_real32_r0, mean_and_stddev_real32_r1, &
    & mean_and_stddev_real64_r0, mean_and_stddev_real64_r1, &
    & mean_and_stddev_int32_r0,  mean_and_stddev_int32_r1,  &
    & mean_and_stddev_int64_r0,  mean_and_stddev_int64_r1

  procedure, public :: mean_and_standard_deviation_per_level => atlas_functionspace_Nodes__mean_and_stddev_per_level

#ifdef FORTRAN_SUPPORTS_FINAL
  final :: atlas_functionspace_Nodes__final
#endif


END TYPE atlas_functionspace_Nodes

interface atlas_functionspace_Nodes
  module procedure atlas_functionspace_Nodes__cptr
  module procedure atlas_functionspace_Nodes__mesh_halo
end interface

!------------------------------------------------------------------------------

