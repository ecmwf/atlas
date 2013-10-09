
! module for testing creation of a Grid
module read_joanna_module
  use common_module
  use datastruct_module, only: create_mesh, DataStructure_type
  implicit none
contains
  
  subroutine read_joanna(filename,g)
    implicit none
    character(len=*), intent(in) :: filename
    type(DataStructure_type), intent(inout)    :: g

    integer                      :: nnode
    integer                      :: nedge
    integer                      :: nface
    integer                      :: ncoin
    integer                      :: before1
    integer                      :: before2
    integer                      :: edge_cnt
    integer                      :: pole_edge_cnt
    integer                      :: internal_edge_cnt

    integer                      :: idx
    real(kind=jprb)              :: vol, Sx, Sy
    integer                      :: ip1,ip2

    integer :: jedge
    integer :: jnode

    call log_info( "Reading mesh "//filename )

    open(5,file=filename,access='sequential',status='old')

    read(5,*) nnode, nedge, nface, ncoin, before2, before1  !nbefore1<nbefore2

    call create_mesh(nnode,nedge,g)

    do jnode = 1, g%nb_nodes
      read(5,*) g%coordinates(jnode,1), g%coordinates(jnode,2), &
        & g%dual_volumes(jnode)
    end do
    g%internal_mesh%nodes_coordinates = g%coordinates

    g%nb_internal_edges = nedge - (before2-before1+1)
    g%nb_pole_edges = before2-before1+1
    g%nb_ghost_nodes = ncoin

    edge_cnt = 0
    pole_edge_cnt = 0
    internal_edge_cnt = 0
    allocate(g%pole_edges(g%nb_pole_edges))
    allocate(g%internal_edges(g%nb_internal_edges))

    do jedge = 1, g%nb_edges
      read(5,*) idx, ip1, ip2
      g%edges(jedge,1) = ip1
      g%edges(jedge,2) = ip2

      if (jedge<before1 .or. jedge>before2) then
        internal_edge_cnt = internal_edge_cnt + 1
        g%internal_edges(internal_edge_cnt) = jedge
      else
        pole_edge_cnt = pole_edge_cnt + 1
        g%pole_edges(pole_edge_cnt) = jedge
      end if
    end do
    g%internal_mesh%faces = g%edges

    do jedge = 1,g%nb_edges
      read(5,*) idx, Sx, Sy
      g%dual_normals(jedge,:) = [Sx,Sy]
    enddo

    allocate(g%ghost_nodes(g%nb_ghost_nodes,2))
    do jnode = 1,g%nb_ghost_nodes
      read(5,*) idx, ip1, ip2
      g%ghost_nodes(jnode,1) = ip1
      g%ghost_nodes(jnode,2) = ip2
    end do

    close(5)
    call log_info( "Done reading mesh "//filename )
  end subroutine read_joanna


end module read_joanna_module