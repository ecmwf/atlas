! ===================================================================
! shallow_water program
! ---------------------
! This program solves the shallow water equations on a sphere,
! initialised with the Rossby-Haurwitz waves.
! ===================================================================
program shallow_water

  use common_module
  use parallel_module, only: parallel_init, parallel_finalise, nproc, nthread
  use gmsh_module, only: write_gmsh_mesh_2d, write_gmsh_fields
  use grib_module, only: write_grib

  use joanna_module, only: read_joanna, write_results_joanna
  use datastruct_module,  only: DataStructure_type, mark_output
  use shallow_water_module, only: &
    & setup_shallow_water, &
    & set_state_rossby_haurwitz, &
    & set_state_zonal_flow, &
    & set_topography_mountain, &
    & set_time_step, &
    & propagate_state, &
    & EQS_MOMENTUM, EQS_VELOCITY
  use mpdata_module, only: &
    MPDATA_GAUGE, MPDATA_STANDARD
  implicit none

  ! Configuration parameters
  real(kind=jprw), parameter :: dt = 20.              ! solver time-step
  integer, parameter         :: nb_steps = 15         ! Number of propagations
  integer, parameter         :: hours_per_step = 24   ! Propagation time
  integer, parameter         :: order = 2             ! Order of accuracy
  integer, parameter         :: scheme = MPDATA_STANDARD
  integer, parameter         :: eqs_type = EQS_MOMENTUM
  logical, parameter         :: write_itermediate_output = .True.

  ! Declarations
  type(DataStructure_type) :: dstruct
  real(kind=jprw), parameter :: hours = 3600.     ! in seconds
  real(kind=jprw), parameter :: days  = 24.*hours ! in seconds
  integer :: jstep = 0
  type(Timer_type) :: wallclock_timer, step_timer

  ! Execution
  call parallel_init()

  call set_log_level(LOG_LEVEL_INFO)
  call set_log_proc(0)
  
  call log_info("Program shallow_water start")
  call read_joanna("data/meshvol.d","data/rtable_lin_T255.d", dstruct)
  call write_gmsh_mesh_2d(dstruct,"data/mesh.msh")

  call setup_shallow_water(eqs_type,dstruct)

  call set_topography_mountain(16000._jprw,dstruct)

  !call set_state_rossby_haurwitz(dstruct)
  call set_state_zonal_flow(dstruct)

  call set_time_step( dt )

  call log_info( "+------------------------------+" )
  call log_info( "| Simulation summary           |" )
  call log_info( "+-------------------+----------+" )
  call log_info( "| MPI Tasks         | "//trim(str(nproc,'(I8)'))//" |" )
  call log_info( "| OMP Threads       | "//trim(str(nthread,'(I8)'))//" |" )
  call log_info( "| glb_nb_nodes      | "//trim(str(dstruct%glb_nb_nodes,'(I8)'))//" |" )
  call log_info( "| glb_nb_edges      | "//trim(str(dstruct%glb_nb_edges,'(I8)'))//" |" )
  call log_info( "| nb_nodes[0]       | "//trim(str(dstruct%nb_nodes,'(I8)'))//" |" )
  call log_info( "| nb_edges[0]       | "//trim(str(dstruct%nb_edges,'(I8)'))//" |" )
  call log_info( "| order             | "//trim(str(order,'(I8)'))//" |" )
  call log_info( "| time step         | "//trim(str(dt,'(F8.1)'))//" |" )
  call log_info( "| total time  (hrs) | "//trim(str(nb_steps*hours_per_step,'(I8)'))//" |" )
  call log_info( "| output rate (hrs) | "//trim(str(hours_per_step,'(I8)'))//" |" )
  call log_info( "+-------------------+----------+ ")

  call mark_output("topography",dstruct)
  call mark_output("height",dstruct)
  call mark_output("depth",dstruct)
  call mark_output("velocity",dstruct)
  call mark_output("momentum",dstruct)
  call mark_output("forcing",dstruct)

  call write_fields()

  call wallclock_timer%start()

  do jstep=1,nb_steps 

    call step_timer%start()
    call propagate_state( hours_per_step*hours, order, scheme, dstruct)
    !call propagate_state( 7740._jprw, order, scheme, dstruct)
    !call propagate_state( dt, order, scheme, dstruct)

    write (log_str, '(A,I3,A,A,F8.2,A,F8.2,A)') &
      & "Propagated to ",jstep*hours_per_step," hours.", &
      & "     step-time = ",step_timer%elapsed(),&
      & "     wall-time = ",wallclock_timer%elapsed(), new_line('A')
    call log_info()
    
    if (write_itermediate_output) call write_fields

  end do ! steps

  ! Write last step anyway if intermediate output is disabled
  if (.not. write_itermediate_output) call write_fields()

  !call write_results_joanna("data/results.d","data/rtable_lin_T255.d",dstruct)

  call log_info("Program shallow_water exit")
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
    call write_grib(dstruct,filename)
    call wallclock_timer%resume()
  end subroutine write_fields

end program shallow_water
