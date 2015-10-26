! (C) Copyright 2013-2015 ECMWF.


TYPE, extends(atlas_object) :: atlas_mesh_Nodes
contains
procedure, public :: size => Nodes__size
procedure, public :: resize => Nodes__resize
procedure, public :: add => Nodes__add
procedure, public :: remove_field => Nodes__remove_field
procedure, private :: field_by_idx  => Nodes__field_by_idx
procedure, public :: field_by_name => Nodes__field_by_name
generic, public :: field => &
    & field_by_idx, &
    & field_by_name
procedure, public :: nb_fields => Nodes__nb_fields
procedure, public :: has_field => Nodes__has_field
procedure, public :: metadata => Nodes__metadata
procedure, public :: str => Nodes__str

procedure, public :: lonlat => Nodes__lonlat
procedure, public :: global_index => Nodes__global_index
procedure, public :: remote_index => Nodes__remote_index
procedure, public :: partition => Nodes__partition
procedure, public :: ghost => Nodes__ghost

procedure, public :: delete => atlas_mesh_Nodes__delete


END TYPE
