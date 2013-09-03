

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
  real                       :: dt
  integer                    :: i
  type(ShallowWaterModel)    :: shallow_water

  ! Execution
  ! ---------
  
  allocate( Grid_class :: grid )

  call read_joanna_mesh(grid,"data/meshvol.d")
  call read_joanna_fields(grid,"data/meshvol.d") 

  call shallow_water%init(grid)
  call shallow_water%set_state_rossby_haurwitz()

  shallow_water%solver%dt_stability = 2
  dt = 3
  do i=1,2
    call shallow_water%solve_time_step( dt )
    write(0,*) "Completed time step. Time: ",shallow_water%state%time
  end do
  
  !call write_gmsh(g)
  call write_gmsh_state(shallow_water%state,"data/results.msh")

  ! Destruction
  ! -----------
  ! Recursively deallocate the grid, functionspaces, fields, ...
  call grid%destruct()

end program main