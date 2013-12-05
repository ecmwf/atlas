
! module for testing creation of a mesh
module gmsh_module
  use parallel_module, only: myproc, nproc, parallel_barrier, Comm_type
  use datastruct_module
  use datastruct
  use common_module, only : jprw,log_info,XX,YY
  implicit none
contains
    
  subroutine write_gmsh_nodal_field(field, comm)
    type(Field_type), intent(in) :: field
    type(Comm_type), intent(in) :: comm
    
    real(kind=jprw), allocatable :: glb_field(:,:)
    integer :: jnode, glb_rows

    glb_rows = comm%glb_field_size()

    call comm%gather( field%data2(), glb_field )

    if (myproc .eq. 0) then
      write(51,'(A)')"$NodeData"
      write(51,*) 1                     ! one string tag:
      write(51,'(A)') '"'//field%name()//'"'      ! the name of the view ("A scalar view")
      write(51,*) 1                     ! one real(kind=jprw) tag:
      write(51,*) 0.0                   ! the time value (0.0)
      write(51,*) 3                     ! three integer tags:
      write(51,*) 0                     ! the time step (0; time steps always start at 0)
      if (field%nb_vars() == 1) write(51,*) 1    ! 1-component (scalar) field
      if (field%nb_vars() == 2) write(51,*) 3    ! 3-component (vector) field
      write(51,*) size(glb_field,1)          ! number of associated nodal values

      if (field%nb_vars() == 1) then
        do jnode=1,glb_rows
          write(51,'(1I8,E18.10)') jnode, glb_field(jnode,1)
        enddo
      endif
      if (field%nb_vars() == 2) then
        do jnode=1,glb_rows
          write(51,'(1I8,E18.10,E18.10,F3.0)') jnode, glb_field(jnode,1), glb_field(jnode,2), 0.0
        enddo
      endif
      write(51,'(A)')"$EndNodeData"
    end if
  end subroutine write_gmsh_nodal_field
  
  subroutine write_gmsh_mesh(dstruct,filename)
    class(DataStructure_type), intent(in) :: dstruct
    character(len=*), intent(in) :: filename
    character(len=1024) :: procfile
    integer :: jedge, jnode
    real(kind=jprw) :: r, phi, theta, x, y, z, pi
    real(kind=jprw), pointer :: coords(:,:)
    pi = acos(-1._jprw)

    coords => vector_field("coordinates",dstruct)

    call log_info( "Writing Gmsh file "//trim(filename) )
    
    call parallel_barrier()

    write(procfile,'(A,A,I2.2)') trim(filename),'.P',myproc
    open(50,file=trim(procfile),access='sequential',status='REPLACE')
    
    write(50,'(A)')"$MeshFormat"
    write(50,'(A)')"2.2 0 8"
    write(50,'(A)')"$EndMeshFormat"
    write(50,'(A)')"$Nodes"
    write(50,*) dstruct%nb_nodes
    do jnode=1, dstruct%nb_nodes
      r     = 6371.
      phi   = coords(jnode,XX)
      theta = -(coords(jnode,YY)+pi/2.)
      x = r*sin(theta)*cos(phi)
      y = r*sin(theta)*sin(phi)
      z = r*cos(theta)

      write(50,'(1I8,F18.10,F18.10,F18.10)')  dstruct%nodes_glb_idx(jnode), x,y,z
      !write(50,'(1I8,F18.10,F18.10,F18.10)')  mesh%nodes%glb_idx(jnode), &
      !  & mesh%nodes_coordinates(jnode,1),mesh%nodes_coordinates(jnode,2), 0.
    enddo
    write(50,'(A)')"$EndNodes"
    write(50,'(A)')"$Elements"
    write(50,*) dstruct%nb_edges
    do jedge=1, dstruct%nb_edges
      ! element-number  type(1=lineP1)  nb_tags(=3)  tag1(=physical-group)  tag2(=elementary-group) tag4(=nb_partitions) tag3(=partition) [nodes]
      write(50,'(1I8,5I2,I3, 2I8)')  dstruct%edges_glb_idx(jedge), 1, 4, 1, 1, 1, dstruct%edges_proc(jedge), &
        & dstruct%nodes_glb_idx( dstruct%edges(jedge,1) ), dstruct%nodes_glb_idx( dstruct%edges(jedge,2) )
    enddo
    write(50,'(A)')"$EndElements"
    close(50)
  end subroutine write_gmsh_mesh

  subroutine write_gmsh_fields(dstruct,filename)
    class(DataStructure_type), intent(in) :: dstruct
    character(len=*), intent(in) :: filename
    type(FieldSet_type) :: fieldset
    type(Field_type) :: field
    type(FunctionSpace_type) :: function_space
    integer :: jfield

    call parallel_barrier()

    call log_info( "Writing Gmsh file "//trim(filename) )
    if (myproc .eq. 0) then
      open(51,file=filename,access='sequential',status='REPLACE')
      write(51,'(A)')"$MeshFormat"
      write(51,'(A)')"2.2 0 8"
      write(51,'(A)')"$EndMeshFormat"
    end if

    fieldset = dstruct%output_fields
    do jfield=1,fieldset%size()
      field = fieldset%field(jfield)
      function_space = field%function_space()
      if (function_space%name() == "nodes") then
        call write_gmsh_nodal_field(field,dstruct%nodes_comm)
      end if
    end do
      
    if (myproc .eq. 0) then
      close(51)
    end if
  end subroutine write_gmsh_fields

end module gmsh_module