

! Main program
program main
  use grid_module
  use read_joana_module
  use gmsh_module
  use shallow_water_module
  implicit none
  
  ! Declarations
  ! ------------
  type(Grid)                        :: g
  class(FunctionSpace),  pointer    :: vertices
  class(FunctionSpace),  pointer    :: faces
  class(Field),  pointer            :: V
  class(Field),  pointer            :: S
  real                              :: dt
  integer                           :: f
  integer                           :: i
  type(ShallowWaterModel) :: shallow_water

  ! Execution
  ! ---------
  
  call read_joanna_mesh(g,"meshvol.d")

  ! this will be replaced by algorithm to compute 
  ! dual volume and dual face normals
  call read_joanna_fields(g,"meshvol.d") 

  vertices => g%function_space("vertices")
  faces    => g%function_space("faces")

  V => vertices%field("dual_volume")
  S => faces%field("dual_face_normal")

  call shallow_water%init(vertices)
  call init_state_rossby_haurwitz(shallow_water%state)

  ! Do 10 time steps of 6 hours in seconds
  shallow_water%solver%dt_stability = 3752
  dt = 6*60*60
  do i=1,10
    call shallow_water%solve_time_step( dt )
    write(0,*) "Completed time step. Time: ",shallow_water%state%time
  end do

  write(0,*) "Shallow Water State Time = ",shallow_water%state%time
  write(0,*) "Shallow Water State Fields:"
  do f=1,size(shallow_water%state%fields)
    write(0,*) " - ",shallow_water%state%fields(f)%ptr%name
  end do
  
  call write_gmsh_state(shallow_water%state,"results.msh")

  ! Destruction
  ! -----------
  ! Recursively deallocate the grid, functionspaces, fields, ...
  call g%destruct()

end program main