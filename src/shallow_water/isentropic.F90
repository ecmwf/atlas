! ===================================================================
! isentropic program
! ---------------------
! This program solves the shallow water equations on a sphere,
! initialised with the Rossby-Haurwitz waves.
! ===================================================================
program isentropic

  use common_module
  use parallel_module, only: parallel_init, parallel_finalise, nproc, nthread
  use gmsh_module, only: write_gmsh_mesh_3d, write_gmsh_fields
  !use grib_module, only: write_grib

  use joanna_module, only: read_joanna_3d, write_results_joanna
  use datastruct_module,  only: DataStructure_type, mark_output
  use isentropic_module, only: &
    & setup_isentropic, &
    & set_state_zonal_flow, &
    & set_topography, &
    & set_time_step, &
    & propagate_state
  implicit none

  ! Configuration parameters
  real(kind=jprw) :: dt = 1.              ! solver time-step
  integer         :: nb_steps = 1         ! Number of propagations
  integer         :: hours_per_step = 1   ! Propagation time
  logical         :: write_itermediate_output = .True.
  integer         :: nb_levels = 5

  ! Declarations
  type(DataStructure_type) :: dstruct
  real(kind=jprw), parameter :: hours = 3600.     ! in seconds
  real(kind=jprw), parameter :: days  = 24.*hours ! in seconds
  integer :: jstep = 0
  type(Timer_type) :: wallclock_timer, step_timer

  ! Execution
  call parallel_init()

  call set_log_level(LOG_LEVEL_DEBUG)
  call set_log_proc(0)
  
  call log_info("Program isentropic start")
  call read_joanna_3d("data/meshvol.d","data/rtable_lin_T255.d", nb_levels, dstruct)

  call log_info( "+------------------------------+" )
  call log_info( "| Simulation summary           |" )
  call log_info( "+-------------------+----------+" )
  call log_info( "| MPI Tasks         | "//trim(str(nproc,'(I8)'))//" |" )
  call log_info( "| OMP Threads       | "//trim(str(nthread,'(I8)'))//" |" )
  call log_info( "| glb_nb_nodes      | "//trim(str(dstruct%glb_nb_nodes,'(I8)'))//" |" )
  call log_info( "| glb_nb_edges      | "//trim(str(dstruct%glb_nb_edges,'(I8)'))//" |" )
  call log_info( "| nb_nodes[0]       | "//trim(str(dstruct%nb_nodes,'(I8)'))//" |" )
  call log_info( "| nb_edges[0]       | "//trim(str(dstruct%nb_edges,'(I8)'))//" |" )
  call log_info( "| nb_levels         | "//trim(str(dstruct%nb_levels,'(I8)'))//" |" )
  call log_info( "| time step         | "//trim(str(dt,'(F8.1)'))//" |" )
  call log_info( "| total time  (hrs) | "//trim(str(nb_steps*hours_per_step,'(I8)'))//" |" )
  call log_info( "| output rate (hrs) | "//trim(str(hours_per_step,'(I8)'))//" |" )
  call log_info( "+-------------------+----------+ ")

  !call write_gmsh_mesh_3d(dstruct,"data/mesh3d.msh")

  call setup_isentropic(dstruct)

  call set_topography(dstruct)

  !call set_state_rossby_haurwitz(dstruct)
  call set_state_zonal_flow(dstruct)

  call set_time_step( dt )

  call mark_output("topography",dstruct)
  call mark_output("height",dstruct)
  call mark_output("depth",dstruct)
  call mark_output("velocity",dstruct)
  call mark_output("pressure",dstruct)
  call mark_output("forcing",dstruct)
  call mark_output("montgomery_potential",dstruct)

  call write_fields()
  call wallclock_timer%start()

  do jstep=1,nb_steps 

    call step_timer%start()
    !call propagate_state( hours_per_step*hours, dstruct)
    call propagate_state( dt, dstruct)

    write (log_str, '(A,I3,A,A,F8.2,A,F8.2,A)') &
      & "Propagated to ",jstep*hours_per_step," hours.", &
      & "     step-time = ",step_timer%elapsed(),&
      & "     wall-time = ",wallclock_timer%elapsed(), new_line('A')
    call log_info()
    
    if (write_itermediate_output) call write_fields

  end do ! steps

  ! Write last step anyway if intermediate output is disabled
  if (.not. write_itermediate_output) call write_fields()

  call log_info("Program isentropic exit")
  call parallel_finalise()


contains

  
  ! Subroutine write_fields
  ! -----------------------
  ! This routine creates a gmsh and a grib file for every step
  ! Don't worry about its implementation, it is likely to change in the future

  subroutine write_fields
    character(len=1024) :: filename
    call wallclock_timer%pause()
    write (filename, "(A,I2.2,A,I2.2)") "data/fields",jstep,".msh"
    call write_gmsh_fields(dstruct,filename)
    write (filename, "(A,I2.2,A)") "data/fields",jstep,".grib"
    !call write_grib(g,filename)
    call wallclock_timer%resume()
  end subroutine write_fields

end program isentropic
