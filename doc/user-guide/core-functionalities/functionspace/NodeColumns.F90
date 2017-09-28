program main
use, intrinsic :: iso_c_binding, only: c_double
use atlas_module
implicit none
integer, parameter                    :: wp = c_double
character(len=1024)                   :: string
character(len=1024)                   :: gridID
character(len=32)                     :: checksum
type(atlas_StructuredGrid)               :: grid
type(atlas_mesh)                      :: mesh
type(atlas_meshgenerator)             :: meshgenerator
type(atlas_Output)                    :: gmsh
type(atlas_functionspace_NodeColumns) :: fs_nodes
type(atlas_mesh_Nodes)                :: meshnodes
type(atlas_Field)                     :: field_scalar1
type(atlas_Field)                     :: field_scalar2
type(atlas_Field)                     :: field_vector1
type(atlas_Field)                     :: field_vector2
type(atlas_Field)                     :: field_tensor1
type(atlas_Field)                     :: field_tensor2
type(atlas_Field)                     :: lonlatField
type(atlas_Field)                     :: field_global
type(atlas_FieldSet)                  :: fields
integer                               :: nb_nodes, jnode
integer                               :: nb_levels = 10
integer                               :: halo_size = 1
type(atlas_Field)                     :: global, scal
real(wp), pointer                     :: scalar1(:)
real(wp), pointer                     :: lonlat(:,:)
real(wp)                              :: minimum, maximum
real(wp)                              :: sum, oisum
real(wp)                              :: mean, stddev
real(wp), allocatable                 :: minimumv(:), maximumv(:)
real(wp), allocatable                 :: sumv(:), oisumv(:)
real(wp), allocatable                 :: meanv(:), stddevv(:)
integer(ATLAS_KIND_GIDX)              :: glb_idx
integer(ATLAS_KIND_GIDX), allocatable :: glb_idxv (:)

! Variables for scalar1 field definition
real(wp), parameter :: rpi = 2._wp * asin(1._wp)
real(wp), parameter :: deg2rad = rpi / 180._wp
real(wp), parameter :: zlatc = 0._wp * rpi
real(wp), parameter :: zlonc = 1._wp * rpi
real(wp), parameter :: zrad  = 2._wp * rpi / 9._wp
real(wp)            :: zdist, zlon, zlat;

call atlas_library%initialise()

! Generate global classic reduced Gaussian grid
gridID = "N32"
grid = atlas_StructuredGrid(gridID)

! Generate mesh associated to structured grid
meshgenerator = atlas_MeshGenerator()
mesh          = meshgenerator%generate(grid)

! Generate functionspace associated to mesh
fs_nodes      = atlas_functionspace_NodeColumns(mesh, halo=halo_size, levels=nb_levels)

! Note on field generation
field_scalar1 = fs_nodes%create_field( &
               & name="scalar1", &
               & kind=atlas_real(wp), &
               & levels=0 )
field_scalar2 = fs_nodes%create_field( &
               & name="scalar2", &
               & kind=atlas_real(wp) )
field_vector1 = fs_nodes%create_field( &
               & name="vector1", &
               & kind=atlas_real(wp), &
               & variables=2 )
field_vector2 = fs_nodes%create_field( &
               & name="vector2", &
               & kind=atlas_real(wp), &
               & variables=2 )
field_tensor1 = fs_nodes%create_field( &
               & name="tensor1", &
               & kind=atlas_real(wp), &
               & levels=0, &
               & variables=2*2 )
field_tensor2 = fs_nodes%create_field( &
               & name="tensor2", &
               & kind=atlas_real(wp), &
               & variables=2*2 )
!........!
! Number of nodes in the mesh
! (different from number of points on a grid!)
meshnodes     = fs_nodes%nodes()
nb_nodes      = fs_nodes%nb_nodes()

! Retrieve lonlat field to calculate scalar1 function
call field_scalar1%data(scalar1)
lonlatField = meshnodes%lonlat()
call lonlatField%data(lonlat)

do jnode=1,nb_nodes
  zlon = lonlat(1,jnode) * deg2rad
  zlat = lonlat(2,jnode) * deg2rad

  zdist = 2._wp * sqrt((cos(zlat) * sin((zlon-zlonc)/2._wp)) *  &
           & (cos(zlat) * sin((zlon-zlonc)/2._wp)) + &
           &  sin((zlat-zlatc)/2._wp) * sin((zlat-zlatc)/2._wp))

  scalar1(jnode) = 0._wp;
  if (zdist < zrad) then
    scalar1(jnode) = 0.5_wp * (1._wp + cos(rpi*zdist/zrad));
  endif
enddo

! Write mesh and field in gmsh format for visualization
gmsh = atlas_output_Gmsh("mesh.msh")
call gmsh%write(mesh)
gmsh = atlas_output_Gmsh("scalar1.msh")
call gmsh%write(field_scalar1)
!........!
! Halo exchange
call fs_nodes%halo_exchange(field_scalar1)

checksum = fs_nodes%checksum(field_scalar1)
write(string, *) checksum
call atlas_log%info(string)

! Create a global field
field_global = fs_nodes%create_field(name="global", template=field_scalar1, global=.true.)

! Gather operation
call fs_nodes%gather(field_scalar1, field_global);

write(string, *) "local nodes          = ", fs_nodes%nb_nodes()
call atlas_log%info(string)

write(string, *) "grid points          = ", grid%size()
call atlas_log%info(string)

write(string, *) "field_global.shape(1) = ", field_global%shape(1)
call atlas_log%info(string)

! Scatter operation
call fs_nodes%scatter(field_global, field_scalar1)

! Halo exchange and checksum
call fs_nodes%halo_exchange(field_scalar1);
checksum = fs_nodes%checksum(field_scalar1);
write(string, *) checksum
call atlas_log%info(string)

! FieldSet checksum
fields = atlas_FieldSet("")
call fields%add(field_scalar1);
call fields%add(field_global);
checksum = fs_nodes%checksum(fields);
write(string, *) checksum
call atlas_log%info(string)
!........!
! Operations

! Minimum and maximum
call fs_nodes%minimum(field_scalar1, minimum)
call fs_nodes%maximum(field_scalar1, maximum)
write(string, *) "min = ",minimum, " max = ", maximum;
call atlas_log%info(string)


! Minimum and maximum + location
call fs_nodes%minimum_and_location(field_scalar1, minimum, glb_idx)
write(string,*) "min = ",minimum, " gidx = ", glb_idx
call atlas_log%info(string)
call fs_nodes%maximum_and_location(field_scalar1, maximum, glb_idx)
write(string,*) "max = ",maximum, " gidx = ", glb_idx
call atlas_log%info(string)

! Summation and order indipedent summation
call fs_nodes%sum(field_scalar1, sum)
call fs_nodes%order_independent_sum(field_scalar1, oisum)
write(string,*) "sum = ", sum, " oisum = ", oisum
call atlas_log%info(string)

! Average over number of nodes
call fs_nodes%mean(field_scalar1, mean)
write(string,*) "mean = ", mean
call atlas_log%info(string)

! Average and standard deviation over number of nodes
call fs_nodes%mean_and_standard_deviation(&
                & field_scalar1, mean, stddev)
write(string,*) "mean = ", mean
call atlas_log%info(string)
write(string,*) "stddev = ", stddev
call atlas_log%info(string)

call grid        %final()
call mesh        %final()
call fs_nodes    %final()
call field_scalar1%final()
call field_scalar2%final()
call field_vector1%final()
call field_vector2%final()
call field_tensor1%final()
call field_tensor2%final()
call field_global %final()
call fields      %final()

call atlas_library%finalise()

end program main
