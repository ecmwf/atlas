program main
use iso_c_binding
use atlas_module
implicit none
integer, parameter                    :: wp = C_DOUBLE
character(len=1024)                   :: gridID
character(len=1024)                   :: checksum
type(atlas_ReducedGrid)               :: reducedGrid
type(atlas_mesh)                      :: mesh
type(atlas_meshgenerator)             :: meshgenerator
type(atlas_functionspace_NodeColumns)       :: fs_nodes
type(atlas_mesh_Nodes)                :: meshnodes
type(atlas_Field)                     :: scalarField1
type(atlas_Field)                     :: scalarField2
type(atlas_Field)                     :: vectorField1
type(atlas_Field)                     :: vectorField2
type(atlas_Field)                     :: tensorField1
type(atlas_Field)                     :: tensorField2
type(atlas_Field)                     :: lonlatField
type(atlas_Field)                     :: globalField
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

call atlas_init()

! Generate global reduced grid
call atlas_resource("--grid", "N32", gridID)
reducedGrid = atlas_ReducedGrid(gridID)

! Generate mesh associated to reduced grid
meshgenerator = atlas_reducedgridmeshgenerator()
mesh          = meshgenerator%generate(reducedGrid)

! Generate functionspace associated to mesh
fs_nodes      = atlas_functionspace_NodeColumns(mesh, halo_size)

! Note on field generation
scalarField1 = fs_nodes%create_field("scalarField1", &
               & atlas_real(wp))
scalarField2 = fs_nodes%create_field("scalarField2", &
               atlas_real(wp), nb_levels)
vectorField1 = fs_nodes%create_field("vectorField1", &
               & atlas_real(wp), [2])
vectorField2 = fs_nodes%create_field("vectorField2", &
               atlas_real(wp), nb_levels, [2])
tensorField1 = fs_nodes%create_field("tensorField1", &
               & atlas_real(wp), [2,2])
tensorField2 = fs_nodes%create_field("tensorField2", &
               atlas_real(wp), nb_levels, [2,2])
!........!
! Number of nodes in the mesh
! (different from number of points on a grid!)
meshnodes     = fs_nodes%nodes()
nb_nodes      = fs_nodes%nb_nodes()

! Retrieve lonlat field to calculate scalar1 function
call scalarField1%data(scalar1)
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
call atlas_write_gmsh      (mesh, "mesh.msh")
call atlas_write_gmsh_field(scalarField1, fs_nodes, "scalar1.msh")
!........!
! Halo exchange
call fs_nodes%halo_exchange(scalarField1)

checksum = fs_nodes%checksum(scalarField1)
if (atlas_mpi_rank() == 0) then
  write(6, *) checksum
endif

! Create a global field
globalField = fs_nodes%create_global_field("global", scalarField1)

! Gather operation
call fs_nodes%gather(scalarField1, globalField);
if (atlas_mpi_rank() == 0) then
  write(6, *) "local nodes          = ", fs_nodes%nb_nodes()
  write(6, *) "grid points          = ", reducedGrid%npts()
  write(6, *) "globalField.shape(1) = ", globalField%shape(1)
endif

! Scatter operation
call fs_nodes%scatter(globalField, scalarField1)

! Halo exchange and checksum
call fs_nodes%halo_exchange(scalarField1);
checksum = fs_nodes%checksum(scalarField1);
if (atlas_mpi_rank() == 0) then
  write(6, *) checksum
endif

! FieldSet checksum
fields = atlas_FieldSet("")
call fields%add(scalarField1);
call fields%add(globalField);
checksum = fs_nodes%checksum(fields);
if (atlas_mpi_rank() == 0) then
  write(6, *) checksum;
endif
!........!
! Operations

! Minimum and maximum
call fs_nodes%minimum(scalarField1, minimum)
call fs_nodes%maximum(scalarField1, maximum)
write(atlas_log%msg, *) "min = ",minimum, " max = ", maximum;
call atlas_log%info()


! Minimum and maximum + location
call fs_nodes%minimum_and_location(scalarField1, minimum, glb_idx)
write(atlas_log%msg,*) "min = ",minimum, " gidx = ", glb_idx
call atlas_log%info()
call fs_nodes%maximum_and_location(scalarField1, maximum, glb_idx)
write(atlas_log%msg,*) "max = ",maximum, " gidx = ", glb_idx
call atlas_log%info()

! Summation and order indipedent summation
call fs_nodes%sum(scalarField1, sum)
call fs_nodes%order_independent_sum(scalarField1, oisum)
write(atlas_log%msg,*) "sum = ", sum, " oisum = ", oisum
call atlas_log%info()

! Average over number of nodes
call fs_nodes%mean(scalarField1, mean)
write(atlas_log%msg,*) "mean = ", mean
call atlas_log%info()

! Average and standard deviation over number of nodes
call fs_nodes%mean_and_standard_deviation(&
                & scalarField1, mean, stddev)
write(atlas_log%msg,*) "mean = ", mean
call atlas_log%info()
write(atlas_log%msg,*) "stddev = ", stddev
call atlas_log%info()

call reducedGrid %final()
call mesh        %final()
call fs_nodes    %final()
call scalarField1%final()
call scalarField2%final()
call vectorField1%final()
call vectorField2%final()
call tensorField1%final()
call tensorField2%final()
call globalField %final()
call fields      %final()

call atlas_finalize()

end program main
