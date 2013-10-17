! ===================================================================
! shallow_water program
! ---------------------
! This program solves the shallow water equations on a sphere,
! initialised with the Rossby-Haurwitz waves.
! ===================================================================
program shallow_water

  use common_module
  use parallel_module
  use gmsh_module, only: write_gmsh_mesh, write_gmsh_state
  use grib_module, only: write_grib

  use read_joanna_module, only: read_joanna
  use datastruct_module,  only: &
    & DataStructure_type, &
    & scalar_field, &
    & vector_field, &
    & create_scalar_field_in_nodes, &
    & create_vector_field_in_nodes, &
    & synchronise
  use shallow_water_module, only: &
    & setup_shallow_water, &
    & set_state_rossby_haurwitz, &
    & set_time_step, &
    & propagate_state
  implicit none

  ! Configuration parameters
  real(kind=jprb) :: dt = 20.              ! solver time-step
  integer         :: nb_steps = 24         ! Number of propagations
  integer         :: hours_per_step = 1   ! Propagation time
  logical         :: write_itermediate_output = .True.

  ! Declarations
  type(DataStructure_type) :: g
  real(kind=jprb), parameter :: hours = 3600.     ! in seconds
  real(kind=jprb), parameter :: days  = 24.*hours ! in seconds
  integer :: jstep = 0
  type(Timer_type) :: wallclock_timer, step_timer

  real(kind=jprb), pointer :: test_field(:,:)
  integer :: jnode

  ! Execution
  call parallel_init()

  call set_log_level(LOG_LEVEL_INFO)
  call set_log_proc(0)
  
  call log_info("Program shallow_water start")
  call read_joanna("data/meshvol.d","data/rtable_lin_T255.d", g)

  call write_gmsh_mesh(g%internal_mesh,"data/mesh.msh")

  call setup_shallow_water(g)
  call set_state_rossby_haurwitz(g)
  call set_time_step( dt )

  call log_info( "+------------------------------+" )
  call log_info( "| Simulation summary           |" )
  call log_info( "+-------------------+----------+" )
  call log_info( "| nb_nodes          | "//trim(str(g%nb_nodes,'(I8)'))//" |" )
  call log_info( "| nb_edges          | "//trim(str(g%nb_edges,'(I8)'))//" |" )
  call log_info( "| time step         | "//trim(str(dt,'(F8.1)'))//" |" )
  call log_info( "| total time  (hrs) | "//trim(str(nb_steps*hours_per_step,'(I8)'))//" |" )
  call log_info( "| output rate (hrs) | "//trim(str(hours_per_step,'(I8)'))//" |" )
  call log_info( "+-------------------+----------+ ")


  call write_fields

  call wallclock_timer%start()

  open(101,file='point.d',access='sequential',status='unknown')


  call create_vector_field_in_nodes( "test_field", g)
  test_field => vector_field("test_field",g)

  test_field(:,1) = myproc
  test_field(:,2) = 0.

  do jnode=1,g%nb_nodes
 
    test_field(jnode,2) = g%internal_mesh%nodes%glb_idx(jnode)

    if (g%internal_mesh%nodes%proc(jnode) .ne. myproc) then
      test_field(jnode,:) = -1
    end if

  end do

  call synchronise("test_field",g)

  do jnode=1,g%nb_nodes
    if (test_field(jnode,1) .ne. g%internal_mesh%nodes%proc(jnode) ) then
      write(0,*) "ERROR: proc of node",g%internal_mesh%nodes%glb_idx(jnode)," is wrong"
    end if
    if (test_field(jnode,2) .ne. g%internal_mesh%nodes%glb_idx(jnode) ) then
      write(0,*) "ERROR: idx of node",g%internal_mesh%nodes%glb_idx(jnode)," is wrong"
    end if
  end do

  do jstep=1,nb_steps 

    call step_timer%start()
    call propagate_state( hours_per_step*hours, g)
  
    write (log_str, '(A,I3,A,A,F8.2,A,F8.2,A)') &
      & "Propagated to ",jstep*hours_per_step," hours.", &
      & "     step-time = ",step_timer%elapsed(),&
      & "     wall-time = ",wallclock_timer%elapsed(), new_line('A')
    call log_info()
    
    if (write_itermediate_output) call write_fields

  end do ! steps

  close(101)

  ! Write last step anyway if intermediate output is disabled
  if (.not. write_itermediate_output) call write_fields

  call write_results_joanna("data/results.d")

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
    write (filename, "(A,I2.2,A,I2.2)") "data/fields",jstep,".msh.P",myproc
    call write_gmsh_state(g%fields,filename)
    write (filename, "(A,I2.2,A)") "data/depth",jstep,".grib"
    call write_grib(g,filename)
    call wallclock_timer%resume()
  end subroutine write_fields


  ! Subroutine write_results_joanna
  ! -------------------------------
  ! This routine generates the results.d file as in 
  ! Joanna Szmelter's original code
  subroutine write_results_joanna(filename)
    character(len=*) :: filename
    integer :: jnode, inode, ip1, ip2, iaaa
    real(kind=jprb), pointer :: D(:), Q(:,:), coords(:,:)
    call wallclock_timer%pause()
    call log_info("Writing output "//filename)
    open(11,file=filename,access='sequential',status='unknown')
    D => scalar_field("depth",g)
    Q => vector_field("momentum",g)
    coords => vector_field("coordinates",g)
    write(11,*)g%fields%time
    write(11,*)0
    write(11,'(A)')'point X Y p'
    do jnode=1,g%nb_nodes

      iaaa=0

      ! special treatment for periodic points specific for Lboro plotting 
      do inode=1,g%nb_ghost_nodes
        ip1=g%ghost_nodes(inode,1)
        ip2=g%ghost_nodes(inode,2)
        if(jnode==ip2)then
          write(11,*) jnode, coords(jnode,XX), coords(jnode,YY), &
            & D(ip1), Q(ip1,XX)/D(ip1), Q(ip1,YY)/D(ip1), D(ip1)
          iaaa=1
        endif
      enddo

      if(iaaa==0)then
        write(11,*) jnode, coords(jnode,XX), coords(jnode,YY), &
          & D(jnode), Q(jnode,XX)/D(jnode), Q(jnode,YY)/D(jnode), D(jnode)
      endif
    enddo
  end subroutine write_results_joanna

end program shallow_water
