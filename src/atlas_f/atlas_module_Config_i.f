! (C) Copyright 2013-2015 ECMWF.

!------------------------------------------------------------------------------

TYPE, extends(object_type) :: atlas_Config

! Purpose :
! -------
!   *Config* : Container of Config, parameters or attributes
!       The Config are seted as key, value pairs

! Methods :
! -------
!   set : set a new property with given key and value
!   set : Modify a property with given key and value
!   get : Return a property value for given key

! Author :
! ------
!   June-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, private :: set_config => atlas_Config__set_config
  procedure, private :: set_config_list => atlas_Config__set_config_list
  procedure, private :: set_logical => atlas_Config__set_logical
  procedure, private :: set_int32 => atlas_Config__set_int32
  procedure, private :: set_real32 => atlas_Config__set_real32
  procedure, private :: set_real64 => atlas_Config__set_real64
  procedure, private :: set_string => atlas_Config__set_string
  procedure, private :: set_array_int32 => atlas_Config__set_array_int32
  procedure, private :: set_array_int64 => atlas_Config__set_array_int64
  procedure, private :: set_array_real32 => atlas_Config__set_array_real32
  procedure, private :: set_array_real64 => atlas_Config__set_array_real64
  procedure :: has => atlas_Config__has
  generic :: set => set_config, set_config_list, set_logical, set_int32, set_real32, set_real64, &
                    set_string, set_array_int32, set_array_int64, set_array_real32, set_array_real64
  procedure, private :: get_config => atlas_Config__get_config
  procedure, private :: get_config_list => atlas_Config__get_config_list
  procedure, private :: get_int32 => atlas_Config__get_int32
  procedure, private :: get_logical => atlas_Config__get_logical
  procedure, private :: get_real32 => atlas_Config__get_real32
  procedure, private :: get_real64 => atlas_Config__get_real64
  procedure, private :: get_string => atlas_Config__get_string
  procedure, private :: get_array_int32 => atlas_Config__get_array_int32
  procedure, private :: get_array_int64 => atlas_Config__get_array_int64
  procedure, private :: get_array_real32 => atlas_Config__get_array_real32
  procedure, private :: get_array_real64 => atlas_Config__get_array_real64
  generic :: get => get_config, get_config_list, get_int32, get_logical, get_real32, get_real64, &
                    get_string, get_array_int32, get_array_int64, get_array_real32, get_array_real64
  procedure :: json => atlas_Config__json

END TYPE atlas_Config

!------------------------------------------------------------------------------

interface atlas_Config
  module procedure atlas_Config__ctor
  module procedure atlas_Config__ctor_from_file
  module procedure atlas_Config__ctor_from_json
end interface

!------------------------------------------------------------------------------
