
! module for testing creation of a Grid
module gmsh_module
  use grid_module
  use model_module
  implicit none
contains
    
  subroutine write_gmsh_nodal_field(g,name)
    class(Grid_class), intent(inout) :: g
    character(len=*) , intent(in) :: name
    class(FunctionSpace_class),  pointer    :: vertices
    class(Field_class),  pointer            :: F
    integer :: inode
    vertices => g%function_space("vertices")
    F => vertices%field(name)
    write(50,'(A)')"$NodeData"
    write(50,*) 1                     ! one string tag:
    write(50,'(A)') '"'//F%name//'"'      ! the name of the view ("A scalar view")
    write(50,*) 1                     ! one real tag:
    write(50,*) 0.0                   ! the time value (0.0)
    write(50,*) 3                     ! three integer tags:
    write(50,*) 0                     ! the time step (0; time steps always start at 0)
    if (F%cols == 1) write(50,*) 1    ! 1-component (scalar) field
    if (F%cols == 2) write(50,*) 3    ! 3-component (vector) field
    write(50,*) g%nb_nodes            ! number of associated nodal values
    do inode=1,g%nb_nodes
      if (F%cols == 1) write(50,*) inode, F%array(inode,1)
      if (F%cols == 2) write(50,*) inode, F%array(inode,1), F%array(inode,2), 0     
    enddo
    write(50,'(A)')"$EndNodeData"
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
    write(50,*) 1                     ! one real tag:
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

  subroutine write_gmsh(g)
    class(Grid_class), intent(inout) :: g
    integer :: iface, inode
    write(0,*) "Writing Gmsh file meshvol.msh"
    open(50,file='meshvol.msh',access='sequential',status='REPLACE')
    write(50,'(A)')"$MeshFormat"
    write(50,'(A)')"2.2 0 8"
    write(50,'(A)')"$EndMeshFormat"
    write(50,'(A)')"$Nodes"
    write(50,*)g%nb_nodes
    do inode=1,g%nb_nodes
      write(50,*) inode, g%nodes(inode,1), g%nodes(inode,2), 0
    enddo
    write(50,'(A)')"$EndNodes"
    write(50,'(A)')"$Elements"
    write(50,*) g%nb_faces
    do iface=1, g%nb_faces
      ! element-number  type(1=lineP1)  nb_tags(=2)  tag1(=physical-group)  tag2(=elementary-group)  [nodes]
      write(50,*)  iface, 1, 2, 1, 1, g%faces(iface,1), g%faces(iface,2)
    enddo
    write(50,'(A)')"$EndElements"
    call write_gmsh_nodal_field(g,"dual_volume")
    call write_gmsh_face_field(g,"dual_face_normal")
    call write_gmsh_nodal_field(g,"depth")
    call write_gmsh_nodal_field(g,"momentum")
    call write_gmsh_nodal_field(g,"rhs_forcing")
    call write_gmsh_nodal_field(g,"advective_velocity")

    close(50)
  end subroutine write_gmsh

  
  subroutine write_gmsh_state(state_,filename)
    class(State_class), intent(in) :: state_
    character(len=*), intent(in) :: filename
    integer :: iface, inode, ifield
    class(Grid_class), pointer :: grid
    grid => state_%function_space%grid
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
      write(50,*)  iface, 1, 2, 1, 1, grid%faces(iface,1), grid%faces(iface,2)
    enddo
    write(50,'(A)')"$EndElements"
    do ifield=1,size(state_%fields)
      call write_gmsh_nodal_field(grid,state_%fields(ifield)%ptr%name)
    end do    
    close(50)
  end subroutine write_gmsh_state

end module gmsh_module