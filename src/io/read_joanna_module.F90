
! module for testing creation of a Grid
module read_joanna_module
  use common_module
  use datastruct_module, only: create_mesh, Geometry_class
  implicit none
contains
  
  subroutine read_joanna(filename,g)
    implicit none
    character(len=*), intent(in) :: filename
    type(Geometry_class), intent(inout)    :: g

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
    integer                      :: p1,p2

    integer :: iedge
    integer :: inode

    open(5,file=filename,access='sequential',status='old')

    read(5,*) nnode, nedge, nface, ncoin, before2, before1  !nbefore1<nbefore2

    call create_mesh(nnode,nedge,g)

    do inode = 1, g%nb_nodes
      read(5,*) g%coordinates(inode,1), g%coordinates(inode,2), &
        & g%dual_volumes(inode)
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

    do iedge = 1, g%nb_edges
      read(5,*) idx, p1, p2
      g%edges(iedge,1) = p1
      g%edges(iedge,2) = p2

      if (iedge<before1 .or. iedge>before2) then
        internal_edge_cnt = internal_edge_cnt + 1
        g%internal_edges(internal_edge_cnt) = iedge
      else
        pole_edge_cnt = pole_edge_cnt + 1
        g%pole_edges(pole_edge_cnt) = iedge
      end if
    end do

    do iedge = 1,g%nb_edges
      read(5,*) idx, Sx, Sy
      g%dual_normals(iedge,:) = [Sx,Sy]
    enddo

    allocate(g%ghost_nodes(g%nb_ghost_nodes,2))
    do inode = 1,g%nb_ghost_nodes
      read(5,*) idx, p1, p2
      g%ghost_nodes(inode,1) = p1
      g%ghost_nodes(inode,2) = p2
    end do

    close(5)
    write(0,*) "Finished reading mesh"
  end subroutine read_joanna


end module read_joanna_module