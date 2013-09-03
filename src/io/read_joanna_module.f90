
! module for testing creation of a Grid
module read_joana_module
  use grid_module
  implicit none
contains
    
  ! read_joanna_mesh(grid)
  ! -----------------
 
  subroutine read_joanna_mesh(g,filename)
    class(Grid_class), pointer, intent(inout)    :: g
    character(len=*), intent(in) :: filename
    integer                      :: nnode
    integer                      :: nedge
    integer                      :: nface
    integer                      :: ncoin
    integer                      :: before1
    integer                      :: before2
    integer                      :: edge_cnt
    integer                      :: idx
    real                         :: vol
    integer                      :: p1,p2

    integer :: iface
    integer :: inode

    open(5,file=filename,access='sequential',status='old')

    read(5,*) nnode, nedge, nface, ncoin, before2, before1  !nbefore1<nbefore2

    ! Specify element type
    call g%init("LagrangeP1_Triag2D")
    
    g%nb_nodes = nnode
    g%nb_elems = 0 ! This mesh only contains edges *sad*
    g%nb_faces = nedge ! -(before2-before1+1)
    
    ! Fill in nodes and cells
    allocate(g%nodes(g%nb_nodes,g%dimension))
    do inode = 1, g%nb_nodes
      read(5,*) g%nodes(inode,1), g%nodes(inode,2), vol
    enddo

    edge_cnt = 0
    allocate(g%faces(g%nb_faces,2)) ! this 2 should be nb_nodes_per_face
    do iface = 1, g%nb_faces
      read(5,*) idx, p1, p2
      if (.true.) then
      !if (iface<before1 .or. iface>before2) then
        edge_cnt = edge_cnt + 1
        g%faces(edge_cnt,1) = p1
        g%faces(edge_cnt,2) = p2
      endif
    enddo
    g%nb_faces = edge_cnt
    
    close(5)
    
  end subroutine read_joanna_mesh


  subroutine read_joanna_fields(g,filename)
    class(Grid_class), pointer, intent(inout) :: g
    character(len=*), intent(in) :: filename

    class(FunctionSpace_class),  pointer    :: vertices
    class(FunctionSpace_class),  pointer    :: faces

    class(Field_class),           pointer    :: V
    class(Field_class),           pointer    :: S

    integer                   :: nnode
    integer                   :: nedge
    integer                   :: nface
    integer                   :: ncoin
    integer                   :: before1
    integer                   :: before2
    integer                   :: edge_cnt
    integer                   :: idx
    real                      :: vol
    integer                   :: p1,p2
    integer :: dummy_int
    real    :: dummy_real
    integer :: iface
    integer :: inode
    integer :: n
    real :: Sx, Sy
    real :: pi,r,y,sm,hx,hy


    open(5,file=filename,access='sequential',status='old')
    read(5,*) nnode, nedge, nface, ncoin, before2, before1  !nbefore1<nbefore2

    ! Create a scalar field for the dual volume in the vertices
    vertices => new_ContinuousFunctionSpace("vertices","LagrangeP1", g)
    V     => vertices%add_scalar_field("dual_volume")

    do inode = 1, g%nb_nodes
      read(5,*) dummy_real, dummy_real, V%array(inode,1)
    enddo

    do iface = 1, nedge
      read(5,*) dummy_int, dummy_int, dummy_int
    enddo

    ! Setup a function space in the face centres for the face_normal
    faces => new_FaceFunctionSpace("faces", "LagrangeP0", g)
    S     => faces%add_vector_field("dual_face_normal")

    edge_cnt = 0
    do iface = 1,nedge
      read(5,*) dummy_int, Sx, Sy
      if (.true.) then
      !if (iface<before1 .or. iface>before2) then
        edge_cnt = edge_cnt + 1
        S%array(edge_cnt,1) = Sx
        S%array(edge_cnt,2) = Sy
      endif
    enddo

    close(5)

  end subroutine read_joanna_fields

end module read_joana_module