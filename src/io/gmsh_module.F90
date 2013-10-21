
! module for testing creation of a Grid
module gmsh_module
  use parallel_module, only: myproc, nproc
  use grid_module
  use mpi
  use common_module, only : jprb,log_info
  implicit none
contains
    
  subroutine write_gmsh_nodal_field(field)
    type(Field_class), pointer, intent(in) :: field
    
    real(kind=jprb), allocatable :: send_field(:)
    integer, allocatable :: send_glb_idx(:)

    real(kind=jprb), allocatable :: glb_field(:,:)
    integer, allocatable :: glb_idx(:)

    integer :: jnode, ierr, owned_nodes_count, jproc, jdim, idx, recv_total_count
    integer :: recv_counts(0:nproc-1)
    integer :: recv_displs(0:nproc-1)
    integer :: root = 0


    owned_nodes_count = 0
    do jnode=1,field%function_space%nb_nodes
      if (field%function_space%proc(jnode) .eq. myproc) then
        owned_nodes_count = owned_nodes_count+1
      end if
    end do

    allocate( send_glb_idx( owned_nodes_count) )
    allocate( send_field( owned_nodes_count) )

    recv_displs(:) = 0
    recv_counts(:) = 0
    call MPI_GATHER( owned_nodes_count, 1, MPI_INTEGER, recv_counts, 1, MPI_INT, root, MPI_COMM_WORLD, ierr ); 
    if (myproc .eq. root) then
      do jproc=1,nproc-1
        recv_displs(jproc) = recv_displs(jproc-1) + recv_counts(jproc-1)
      end do
    end if
    recv_total_count = sum(recv_counts)
    allocate( glb_idx(recv_total_count) )
    allocate( glb_field(recv_total_count,field%cols) )

    idx = 0
    do jnode=1,field%function_space%nb_nodes
      if (field%function_space%proc(jnode) .eq. myproc) then
        idx = idx + 1
        send_glb_idx(idx) = field%function_space%glb_idx(jnode)
      end if
    end do

    call MPI_GATHERV( send_glb_idx, owned_nodes_count, MPI_INTEGER, &
                    & glb_idx, recv_counts, recv_displs, MPI_INTEGER, &
                    & root, MPI_COMM_WORLD, ierr )

    do jdim=1,field%cols
      idx = 0
      do jnode=1,field%function_space%nb_nodes
        if (field%function_space%proc(jnode) .eq. myproc) then
          idx = idx + 1
          send_field(idx) = field%array(jnode,jdim)
        end if
      end do

      call MPI_GATHERV( send_field, owned_nodes_count, MPI_DOUBLE_PRECISION, &
                      & glb_field(:,jdim), recv_counts, recv_displs, MPI_DOUBLE_PRECISION, &
                      & root, MPI_COMM_WORLD, ierr )
    end do

    if (myproc .eq. 0) then
      write(51,'(A)')"$NodeData"
      write(51,*) 1                     ! one string tag:
      write(51,'(A)') '"'//field%name//'"'      ! the name of the view ("A scalar view")
      write(51,*) 1                     ! one real(kind=jprb) tag:
      write(51,*) 0.0                   ! the time value (0.0)
      write(51,*) 3                     ! three integer tags:
      write(51,*) 0                     ! the time step (0; time steps always start at 0)
      if (field%cols == 1) write(51,*) 1    ! 1-component (scalar) field
      if (field%cols == 2) write(51,*) 3    ! 3-component (vector) field
      write(51,*) recv_total_count          ! number of associated nodal values

      if (field%cols == 1) then
        do jnode=1,recv_total_count
          write(51,'(1I8,E18.10)') glb_idx(jnode), glb_field(jnode,1)
        enddo
      endif
      if (field%cols == 2) then
        do jnode=1,recv_total_count
          write(51,'(1I8,E18.10,E18.10,F3.0)') glb_idx(jnode), glb_field(jnode,1), glb_field(jnode,2), 0.0
        enddo
      endif
      write(51,'(A)')"$EndNodeData"
    end if
  end subroutine write_gmsh_nodal_field

  subroutine write_gmsh_face_field(g,name)
    class(Grid_class), intent(inout) :: g
    character(len=*) , intent(in) :: name
    class(FunctionSpace_class),  pointer    :: faces
    class(Field_class),  pointer            :: F
    integer :: jface
    faces => g%function_space("faces")
    F => faces%field(name)
    write(50,'(A)')"$ElementData"
    write(50,*) 1                     ! one string tag:
    write(50,'(A)') '"'//F%name//'"'      ! the name of the view ("A scalar view")
    write(50,*) 1                     ! one real(kind=jprb) tag:
    write(50,*) 0.0                   ! the time value (0.0)
    write(50,*) 3                     ! three integer tags:
    write(50,*) 0                     ! the time step (0; time steps always start at 0)
    if (F%cols == 1) write(50,*) 1    ! 1-component (scalar) field
    if (F%cols == 2) write(50,*) 3    ! 3-component (vector) field
    write(50,*) g%nb_faces            ! number of associated nodal values
    do jface=1,g%nb_faces
      if (F%cols == 1) write(50,*) jface, F%array(jface,1)
      if (F%cols == 2) write(50,*) jface, F%array(jface,1), F%array(jface,2), 0.0     
    enddo
    write(50,'(A)')"$EndElementData"
  end subroutine write_gmsh_face_field

  
  subroutine write_gmsh_mesh(grid,filename)
    class(Grid_class), intent(inout) :: grid
    character(len=*), intent(in) :: filename
    character(len=1024) :: procfile
    integer :: jface, jnode, jproc, ierr
    real(kind=jprb) :: r, phi, theta, x, y, z, pi
    pi = acos(-1._jprb)
    call log_info( "Writing Gmsh file "//trim(filename) )
    !open(40,file=filename,access='sequential',status='REPLACE')
    !do jproc=0,nproc-1
    !  write(40,'(A,A,A,I2.2,A)') 'Merge "',trim(filename),'.P',jproc,'";'
    !end do
    !close(40)
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    write(procfile,'(A,A,I2.2)') trim(filename),'.P',myproc
    open(50,file=trim(procfile),access='sequential',status='REPLACE')

    write(50,'(A)')"$MeshFormat"
    write(50,'(A)')"2.2 0 8"
    write(50,'(A)')"$EndMeshFormat"
    write(50,'(A)')"$Nodes"
    write(50,*)grid%nb_nodes
    do jnode=1,grid%nb_nodes
      r     = 6371.
      phi   = grid%nodes_coordinates(jnode,1)
      theta = -(grid%nodes_coordinates(jnode,2)+pi/2.)
      x = r*sin(theta)*cos(phi)
      y = r*sin(theta)*sin(phi)
      z = r*cos(theta)

      !write(50,'(1I8,F18.10,F18.10,F18.10)')  grid%nodes%glb_idx(jnode), x,y,z
      write(50,'(1I8,F18.10,F18.10,F18.10)')  grid%nodes%glb_idx(jnode), &
        & grid%nodes_coordinates(jnode,1),grid%nodes_coordinates(jnode,2), 0.
    enddo
    write(50,'(A)')"$EndNodes"
    write(50,'(A)')"$Elements"
    write(50,*) grid%nb_faces
    do jface=1, grid%nb_faces
      ! element-number  type(1=lineP1)  nb_tags(=3)  tag1(=physical-group)  tag2(=elementary-group) tag4(=nb_partitions) tag3(=partition) [nodes]
      write(50,'(1I8,5I2,I3, 2I8)')  grid%faces_glb_idx(jface), 1, 4, 1, 1, 1, grid%faces_proc(jface), &
        & grid%nodes%glb_idx( grid%faces(jface,1) ), grid%nodes%glb_idx( grid%faces(jface,2) )
    enddo
    write(50,'(A)')"$EndElements"
    close(50)
  end subroutine write_gmsh_mesh

  subroutine write_gmsh_state(state,filename)
    class(State_class), intent(in), target :: state
    character(len=*), intent(in) :: filename
    integer :: jface, jnode, ifield, ierr
    type(Field_class), pointer :: field

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    field => state%fields(1)%ptr
    call log_info( "Writing Gmsh file "//trim(filename) )
    if (myproc .eq. 0) then
      open(51,file=filename,access='sequential',status='REPLACE')
      write(51,'(A)')"$MeshFormat"
      write(51,'(A)')"2.2 0 8"
      write(51,'(A)')"$EndMeshFormat"
    end if
    field => state%field("depth")
    call write_gmsh_nodal_field(field)

    field => state%field("momentum")
    call write_gmsh_nodal_field(field)

    !do ifield=1,state%nb_fields
    !  field => state%fields(ifield)%ptr
    !  if (field%function_space%name == "nodes") then
    !    call write_gmsh_nodal_field(field)
    !  endif
    !end do 
      
    if (myproc .eq. 0) then
      close(51)
    end if
  end subroutine write_gmsh_state

end module gmsh_module