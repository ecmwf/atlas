
! module for testing creation of a Grid
module gmsh_module
  use grid_module
  use model_module
  implicit none
contains
    
  subroutine write_gmsh_nodal_field(g,field)
    class(Grid_class), intent(inout) :: g
    class(Field_class),  pointer, intent(in) :: field
    integer :: inode
    write(51,'(A)')"$NodeData"
    write(51,*) 1                     ! one string tag:
    write(51,'(A)') '"'//field%name//'"'      ! the name of the view ("A scalar view")
    write(51,*) 1                     ! one real(kind=jprb) tag:
    write(51,*) 0.0                   ! the time value (0.0)
    write(51,*) 3                     ! three integer tags:
    write(51,*) 0                     ! the time step (0; time steps always start at 0)
    if (field%cols == 1) write(51,*) 1    ! 1-component (scalar) field
    if (field%cols == 2) write(51,*) 3    ! 3-component (vector) field
    write(51,*) g%nb_nodes            ! number of associated nodal values
    do inode=1,g%nb_nodes
      if (field%cols == 1) write(51,*) inode, field%array(inode,1)
      if (field%cols == 2) write(51,*) inode, field%array(inode,1), field%array(inode,2), 0     
    enddo
    write(51,'(A)')"$EndNodeData"
  end subroutine write_gmsh_nodal_field

  subroutine write_gmsh_face_field(g,name)
    class(Grid_class), intent(inout) :: g
    character(len=*) , intent(in) :: name
    class(FunctionSpace_class),  pointer    :: faces
    class(Field_class),  pointer            :: F
    integer :: iface
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
    do iface=1,g%nb_faces
      if (F%cols == 1) write(50,*) iface, F%array(iface,1)
      if (F%cols == 2) write(50,*) iface, F%array(iface,1), F%array(iface,2), 0     
    enddo
    write(50,'(A)')"$EndElementData"
  end subroutine write_gmsh_face_field

  
  subroutine write_gmsh_mesh(grid,filename)
    class(Grid_class), intent(inout) :: grid
    character(len=*), intent(in) :: filename
    integer :: iface, inode
    write(0,*) "Writing Gmsh file ",filename
    open(50,file=filename,access='sequential',status='REPLACE')

    write(50,'(A)')"$MeshFormat"
    write(50,'(A)')"2.2 0 8"
    write(50,'(A)')"$EndMeshFormat"
    write(50,'(A)')"$Nodes"
    write(50,*)grid%nb_nodes
    do inode=1,grid%nb_nodes
      write(50,*) inode, grid%nodes(inode,1), grid%nodes(inode,2), 0
    enddo
    write(50,'(A)')"$EndNodes"
    write(50,'(A)')"$Elements"
    write(50,*) grid%nb_faces
    do iface=1, grid%nb_faces
      ! element-number  type(1=lineP1)  nb_tags(=2)  tag1(=physical-group)  tag2(=elementary-group)  [nodes]
      write(50,'(1I8,4I2,2I8)')  iface, 1, 2, 1, 1, grid%faces(iface,1), grid%faces(iface,2)
    enddo
    write(50,'(A)')"$EndElements"
    close(50)
  end subroutine write_gmsh_mesh

  subroutine write_gmsh_state(state_,filename)
    class(State_class), intent(in) :: state_
    character(len=*), intent(in) :: filename
    integer :: iface, inode, ifield
    class(Grid_class), pointer :: grid
    grid => state_%function_space%grid
    write(0,*) "Writing Gmsh file ",filename
    open(51,file=filename,access='sequential',status='REPLACE')
    write(51,'(A)')"$MeshFormat"
    write(51,'(A)')"2.2 0 8"
    write(51,'(A)')"$EndMeshFormat"
    do ifield=1,size(state_%fields)
      call write_gmsh_nodal_field(grid,state_%fields(ifield)%ptr)
    end do    
    close(51)

  end subroutine write_gmsh_state

end module gmsh_module