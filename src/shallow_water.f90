

! Main program
program main
  use grid_module
  use read_joana_module
  use gmsh_module
  use shallow_water_module
  implicit none
  
  ! Declarations
  ! ------------
  class(Grid_class), pointer :: grid
  real(kind=jprb)            :: dt
  integer                    :: i
  type(ShallowWaterModel)    :: shallow_water
  character(len=1024)        :: filename

  ! Execution
  ! ---------
  
  allocate( Grid_class :: grid )

  call read_joanna_mesh(grid,"data/meshvol.d")
  call read_joanna_fields(grid,"data/meshvol.d") 

  call write_gmsh_mesh(grid,"data/mesh.msh")

  call shallow_water%init(grid)
  call shallow_water%set_state_rossby_haurwitz()

  shallow_water%solver%dt_stability = 50.
 
  dt = 60*60*24 / 24 ! = one hour

  do i=1,24 ! every cycle is one dt
    call shallow_water%solve_time_step( dt )
    write(0,*) "Completed time step. Time: ",shallow_water%state%time
    write (filename, "(A12,I5.5,A4)") "data/results",i,".msh"
    call write_gmsh_state(shallow_water%state,filename)
  end do
  
  ! Destruction
  ! -----------
  ! Recursively deallocate the grid, functionspaces, fields, ...
  call grid%destruct()

end program main