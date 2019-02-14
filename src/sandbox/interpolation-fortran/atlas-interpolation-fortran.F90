! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

program atlas_interpolation_fortran
use atlas_module
implicit none

type(atlas_Grid)                      :: grid_A
type(atlas_Grid)                      :: grid_B
type(atlas_Mesh)                      :: mesh_A
type(atlas_Mesh)                      :: mesh_B
type(atlas_functionspace_NodeColumns) :: fs_A
type(atlas_functionspace_NodeColumns) :: fs_B
type(atlas_Field)                     :: field_A, field_A2
type(atlas_Field)                     :: field_B
type(atlas_MeshGenerator)             :: meshgenerator
type(atlas_Partitioner)               :: partitioner_B
type(atlas_GridDistribution)          :: distribution_B
type(atlas_Output)                    :: gmsh
type(atlas_Config)                    :: interpolation_config
type(atlas_Interpolation)             :: interpolation_AB
type(atlas_Interpolation)             :: interpolation_BA
type(atlas_Trace)                     :: trace

call atlas_library%initialise()

trace = atlas_Trace( "atlas-interpolation-fortran.F90", __LINE__, "Complete execution" )
! Setup a meshgenerator
meshgenerator = atlas_MeshGenerator()

! Generate source mesh
grid_A = atlas_Grid("O32")
mesh_A = meshgenerator%generate(grid_A)

! Generate target mesh based on domain decomposition of source mesh
grid_B = atlas_Grid("N16")
partitioner_B  = atlas_MatchingMeshPartitioner(mesh_A)
distribution_B = partitioner_B%partition(grid_B)
mesh_B = meshgenerator%generate(grid_B,distribution_B)

! Create function spaces for each mesh
fs_A = atlas_functionspace_NodeColumns(mesh_A,halo=1)
fs_B = atlas_functionspace_NodeColumns(mesh_B,halo=2)

! Setup interpolators
interpolation_config = atlas_Config()
call interpolation_config%set("type","finite-element")
interpolation_AB = atlas_Interpolation(interpolation_config,fs_A,fs_B)
interpolation_BA = atlas_Interpolation(interpolation_config,fs_B,fs_A)

! Create fields and initialise source field
field_A  = fs_A%create_field(name="A",  kind=atlas_real(atlas_kind_real64))
field_B  = fs_B%create_field(name="B",  kind=atlas_real(atlas_kind_real64))
field_A2 = fs_A%create_field(name="A2", kind=atlas_real(atlas_kind_real64))
call initialise_field_hill(fs_A, field_A)





! Interpolate from source to target
call interpolation_AB%execute(field_A, field_B)
call fs_B%halo_exchange(field_B)

! Interpolate from target back to source
call interpolation_BA%execute(field_B, field_A2)
call fs_A%halo_exchange(field_A2)





! Output target
gmsh = atlas_output_Gmsh("B.msh",coordinates="xyz")
call gmsh%write(mesh_B)
call gmsh%write(field_B)

! Output source
gmsh = atlas_output_Gmsh("A.msh",coordinates="xyz")
call gmsh%write(mesh_A)
call gmsh%write(field_A)
call gmsh%write(field_A2)

! cleanup
call interpolation_config%final()
call interpolation_AB%final()
call interpolation_BA%final()
call gmsh%final()
call partitioner_B%final()
call distribution_B%final()
call field_A%final()
call field_B%final()
call meshgenerator%final()
call fs_A%final()
call fs_B%final()
call mesh_A%final()
call mesh_B%final()
call grid_A%final()
call grid_B%final()

call trace%final()

call atlas_library%finalise()
contains

subroutine initialise_field_hill(funcspace,field)
  type(atlas_functionspace_NodeColumns), intent(in) :: funcspace
  type(atlas_Field), intent(inout) :: field
  real(atlas_kind_real64), parameter :: M_PI = 3.14159265358979323846
  real(atlas_kind_real64), parameter :: deg2rad = M_PI/180._8
  type(atlas_mesh_Nodes) :: nodes
  type(atlas_Field) :: field_lonlat
  real(atlas_kind_real64), pointer :: value(:), lonlat(:,:)
  integer :: jnode, nb_nodes
  real(atlas_kind_real64) :: lon, lat, c2, c_lon, c_lat, c_rad, dist, s1, s2
  c_lat = 0. * M_PI
  c_lon = 1. * M_PI
  c_rad = 2. * M_PI / 9.
  nodes = funcspace%nodes()
  field_lonlat = nodes%lonlat()
  call field_lonlat%data(lonlat)
  call field%data(value)
  nb_nodes = nodes%size()
  do jnode=1,nb_nodes
    lon = deg2rad * lonlat(1,jnode)
    lat = deg2rad * lonlat(2,jnode)
    c2  = cos(lat)
    s1  = sin( (lon-c_lon)/2. )
    s2  = sin( (lat-c_lat)/2. )
    dist = 2. * sqrt( c2*s1*c2*s1 + s2*s2 )
    if( dist < c_rad ) then
      value(jnode) = 1. + cos(M_PI*dist/c_rad)
    else
      value(jnode) = 0
    endif
  enddo
end subroutine

end program
