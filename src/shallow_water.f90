

! Main program
program main
  use grid_module
  use read_joana_module
  use gmsh_module
  use grib_module
  use shallow_water_module
  implicit none
  
  ! Declarations
  ! ------------
  class(Grid_class), pointer :: grid
  real(kind=jprb)            :: step
  integer :: nb_steps
  integer                    :: istep
  class(ShallowWaterModel), pointer :: shallow_water
  character(len=1024)        :: filename

  ! Execution
  ! ---------
  
  allocate( Grid_class :: grid )

  call read_joanna_mesh(grid,"data/meshvol.d")
  call read_joanna_fields(grid,"data/meshvol.d") 

  call write_gmsh_mesh(grid,"data/mesh.msh")

  allocate( shallow_water )
  call shallow_water%init(grid)
  call shallow_water%set_state_rossby_haurwitz()


  write (filename, "(A,I5.5,A)") "data/results",0,".msh"
  call write_gmsh_state(shallow_water%state,filename)

  write (filename, "(A,I5.5,A)") "data/results",0,".grib"
  call write_grib_state(shallow_water%state,filename)


  shallow_water%solver%dt_stability = 20.
 
  step = 60*60*24 ! = one day
  nb_steps = 15
  do istep=1,nb_steps ! every cycle is one dt

    call shallow_water%solve_time_step( step )

    write(0,*) "Completed time step. Time: ",shallow_water%state%time
    
    write (filename, "(A,I5.5,A)") "data/results",istep,".msh"
    call write_gmsh_state(shallow_water%state,filename)

    write (filename, "(A,I5.5,A)") "data/results",istep,".grib"
    call write_grib_state(shallow_water%state,filename)
  end do
  
  ! Destruction
  ! -----------
  ! Recursively deallocate the grid, functionspaces, fields, ...
  call grid%destruct()

end program main