

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
  character(len=1024)        :: filename
  class(FunctionSpace_class), pointer :: faces
  class(Field_class), pointer :: face_normal

  class(State_class), pointer :: mesh_state


  ! Execution
  ! ---------
  
  allocate( Grid_class :: grid )

  call read_joanna_mesh(grid,"data/meshvol.d")
  call read_joanna_fields(grid,"data/meshvol.d") 

  faces => grid%function_space("faces")
  face_normal => faces%field("dual_face_normal")
  allocate( mesh_state )
  call mesh_state%init("mesh_state",faces)
  call mesh_state%add_field( face_normal )

  call write_gmsh_mesh(grid,"data/mesh.msh")
  call write_gmsh_state(mesh_state,"data/faces.msh")
  call shallow_water%init(grid)
  call shallow_water%set_state_rossby_haurwitz()

  shallow_water%solver%dt_stability = 20.
  dt = 20.
  do i=1,1
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