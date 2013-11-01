
! module for testing creation of a Grid
module joanna_module
  use common_module
  use parallel_module
  use datastruct_module, only: create_mesh, DataStructure_type, vector_field, scalar_field
  use split_globe_module, only : split_globe
  use mpi
  implicit none
contains
  
  subroutine read_joanna(filename,rtable,g)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: rtable
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
    integer                      :: owned_node_cnt = 0
    integer                      :: owned_edge_cnt = 0
    integer                      :: glb_nb_nodes = 0
    integer                      :: glb_nb_edges = 0

    integer                      :: idx,ibound
    real(kind=jprw), pointer     :: coords(:,:), vol(:), S(:,:)
    integer                      :: ip1,ip2

    real(kind=jprw) :: dummy_real, x, y, v, sx, sy
    integer :: dummy_int
    
    integer, allocatable :: node_proc_full(:)
    integer, allocatable :: proc(:)
    integer, allocatable :: node_glb_idx_full(:)
    integer, allocatable :: glb_idx(:)
    integer, allocatable :: keep_edge(:)
    integer, allocatable :: keep_node(:)
    integer, allocatable :: node_loc_idx(:)
    integer, allocatable :: edge_loc_idx(:)
    integer, allocatable :: edges_distr(:)

    integer :: ierr
    integer :: jedge, iedge
    integer :: jnode, inode
    integer :: nb_pole_edges

    call log_info( "Reading mesh "//filename )


    ! First pass of the file to see which parts to store in second pass
    open(5,file=filename,access='sequential',status='old')

    read(5,*) nnode, nedge, nface, ncoin, before2, before1  !nbefore1<nbefore2

    allocate( node_proc_full(nnode) )
    allocate( node_glb_idx_full(nnode) )
    allocate( node_loc_idx(nnode) )
    allocate( edge_loc_idx(nedge) )
    allocate( keep_node(nnode) )
    allocate( keep_edge(nedge) )

    !write(log_str,*) "Domain decomposition in ",nproc," parts, using ",rtable; call log_info()
    call split_globe(rtable,nproc,node_proc_full,node_glb_idx_full)

    keep_node(:) = 0
    keep_edge(:) = 0
    node_loc_idx(:) = -1 ! invalidate any loc_idx
    edge_loc_idx(:) = -1 ! invalidate any loc_idx
    do jnode = 1, nnode
      read(5,*) dummy_real, dummy_real, dummy_real
    end do

    nb_pole_edges = 0
    idx=0
    do jedge = 1, nedge
      read(5,*) idx, ip1, ip2
      if ( (node_proc_full(ip1) .eq. myproc) .or. (node_proc_full(ip2) .eq. myproc) ) then
        keep_edge(jedge) = 1
        keep_node(ip1) = 1
        keep_node(ip2) = 1
        if (.not.(jedge<before1 .or. jedge>before2)) then
          nb_pole_edges = nb_pole_edges + 1
        end if
      end if
    end do

    close(5)

    !call log_info(str(myproc)//" should keep nodes:"//str(sum(keep_node)))
    !call log_info(str(myproc)//" should keep edges:"//str(sum(keep_edge)))

    allocate( proc( sum(keep_node) ) ) 
    allocate( glb_idx( sum(keep_node) ) ) 
    idx = 0
    do jnode=1,nnode
      if (keep_node(jnode).eq.1) then
        idx = idx+1
        node_loc_idx(jnode) = idx
        proc(idx) = node_proc_full(jnode)
        glb_idx(idx) = node_glb_idx_full(jnode)
      end if
    end do

    idx = 0
    do jedge=1,nedge
      if (keep_edge(jedge).eq.1) then
        idx = idx+1
        edge_loc_idx(jedge) = idx
      end if
    end do


    ! Create edge-based unstructured mesh
    call create_mesh( sum(keep_node), sum(keep_edge), proc, glb_idx, g )
    
    allocate( g%internal_mesh%faces_proc( g%nb_edges ) )
    allocate( g%internal_mesh%faces_glb_idx( g%nb_edges ) )

    coords => vector_field("coordinates",g)
    vol    => scalar_field("dual_volumes",g)
    S      => vector_field("dual_normals",g)

    open(5,file=filename,access='sequential',status='old')

    read(5,*) nnode, nedge, nface, ncoin, before2, before1  !nbefore1<nbefore2

    owned_node_cnt = 0
    do jnode = 1, nnode
      read(5,*) x, y, v
      if ( keep_node(jnode).eq.1 ) then
        inode = node_loc_idx(jnode)
        coords(inode,:) = [x, y]
        vol(inode) = v
        if( proc(inode) .eq. myproc ) then
          owned_node_cnt = owned_node_cnt + 1
        end if
      end if
    end do

    call MPI_ALLREDUCE( owned_node_cnt, glb_nb_nodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    g%glb_nb_nodes = glb_nb_nodes
    g%internal_mesh%glb_nb_nodes = glb_nb_nodes
    g%internal_mesh%nodes%glb_nb_nodes = glb_nb_nodes
    g%internal_mesh%nodes_coordinates = coords

    g%nb_pole_edges = nb_pole_edges

    g%nb_ghost_nodes = 0!ncoin

    edge_cnt = 0
    pole_edge_cnt = 0
    internal_edge_cnt = 0
    allocate(g%pole_edges(g%nb_pole_edges))

    owned_edge_cnt = 0
    do jedge = 1, nedge
      read(5,*) idx, ip1, ip2
      if ( keep_edge(jedge).eq.1 ) then
        iedge = edge_loc_idx(jedge)
        g%internal_mesh%faces_glb_idx(iedge) = jedge
        g%internal_mesh%faces_proc(iedge) = proc( node_loc_idx(ip1) )
        g%edges(iedge,1) = node_loc_idx(ip1)
        g%edges(iedge,2) = node_loc_idx(ip2)

        if (jedge<before1 .or. jedge>before2) then
          internal_edge_cnt = internal_edge_cnt + 1
        else
          pole_edge_cnt = pole_edge_cnt + 1
          g%pole_edges(pole_edge_cnt) = iedge
        end if

        if( g%internal_mesh%faces_proc(iedge) .eq. myproc ) then
          owned_edge_cnt = owned_edge_cnt + 1
        end if

      end if
    end do
    g%internal_mesh%faces = g%edges

    call MPI_ALLREDUCE( owned_edge_cnt, glb_nb_edges, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    g%glb_nb_edges = glb_nb_edges


    do jedge = 1,nedge
      read(5,*) idx, sx, sy
      if ( keep_edge(jedge).eq.1 ) then
        iedge = edge_loc_idx(jedge)
        S(iedge,:) = [sx, sy]
      end if
    enddo

    allocate(g%ghost_nodes(g%nb_ghost_nodes,2))
    do jnode = 1,g%nb_ghost_nodes
      read(5,*) idx, ip1, ip2
      g%ghost_nodes(jnode,1) = ip1
      g%ghost_nodes(jnode,2) = ip2
    end do

    allocate(g%nb_neighbours(g%nb_nodes))
    g%nb_neighbours(:) = 0
    do jedge = 1,g%nb_edges
      ip1 = g%edges(jedge,1)
      ip2 = g%edges(jedge,2)
      g%nb_neighbours(ip1) = g%nb_neighbours(ip1)+1
      g%nb_neighbours(ip2) = g%nb_neighbours(ip2)+1
    enddo
    g%max_nb_neighbours = maxval(g%nb_neighbours(:))
    ibound = 0
    do jnode = 1,g%nb_nodes
      if( g%nb_neighbours(jnode) == 1) ibound=ibound+1
    enddo
    !do jnode=0,nproc
    !  call mpi_barrier( MPI_COMM_WORLD, ierr)
    !  if(jnode == myproc) &
    !   write(0,*) 'g%max_nb_neighbours  1',myproc,g%max_nb_neighbours, minval(g%nb_neighbours(:)),&
    !   & g%nb_nodes-ibound,ibound,g%nb_edges
    !enddo
    allocate(g%neighbours(g%max_nb_neighbours,g%nb_nodes))
    allocate(g%my_edges(g%max_nb_neighbours,g%nb_nodes))
    allocate(g%sign(g%max_nb_neighbours,g%nb_nodes))
    g%nb_neighbours(:) = 0
    do jedge = 1,g%nb_edges
      ip1 = g%edges(jedge,1)
      ip2 = g%edges(jedge,2)
      g%nb_neighbours(ip1) = g%nb_neighbours(ip1)+1
      g%neighbours(g%nb_neighbours(ip1),ip1) = ip2
      g%my_edges(g%nb_neighbours(ip1),ip1) = jedge
      g%sign(g%nb_neighbours(ip1),ip1) = 1
      g%nb_neighbours(ip2) = g%nb_neighbours(ip2)+1
      g%neighbours(g%nb_neighbours(ip2),ip2) = ip1
      g%my_edges(g%nb_neighbours(ip2),ip2) = jedge
      g%sign(g%nb_neighbours(ip2),ip2) = -1
    enddo
    !if(myproc == 5) then
    !  write(0,*) 'g%nb_nodes ',g%nb_nodes
    !  do jnode = 1,g%nb_nodes
    !    write(0,*) 'neighbours ',jnode,g%nb_neighbours(jnode),g%neighbours(1:g%nb_neighbours(jnode),jnode)
    !  enddo
    !endif
    close(5)

    allocate( edges_distr(nproc) )
    call MPI_GATHER( g%nb_edges, 1, MPI_INTEGER, &
                   & edges_distr, 1, MPI_INTEGER, &
                   & 0, MPI_COMM_WORLD, ierr );

    if( myproc .eq. 0 ) then
      call log_info('           Edges Load balance per MPI task')
      call plot1d(nproc,10,edges_distr)
    end if

  end subroutine read_joanna


  ! Subroutine write_results_joanna
  ! -------------------------------
  ! This routine generates the results.d file as in 
  ! Joanna Szmelter's original code

  subroutine write_results_joanna(filename,rtable,g)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: rtable
    type(DataStructure_type), intent(inout)    :: g

    integer :: jnode, inode, ip1, ip2, iaaa
    integer, allocatable :: plot_idx(:)
    real(kind=jprw), pointer :: loc_D(:), loc_Q(:,:), loc_coords(:,:)
    real(kind=jprw), allocatable :: D(:), Q(:,:), coords(:,:)
    integer(kind=jpim)              :: ndgl     ! Number of lattitudes
    integer(kind=jpim), allocatable :: nloen(:) ! Number of longitude points for each latitude
    integer(kind=jpim) :: jlat,idum, nb_nodes, inode_grwch, jglb

    call log_info("Writing output "//filename)

    loc_D => scalar_field("depth",g)
    loc_Q => vector_field("momentum",g)
    loc_coords => vector_field("coordinates",g)

    call g%internal_mesh%nodes%comm%gather(loc_D,D)
    call g%internal_mesh%nodes%comm%gather(loc_Q,Q)
    call g%internal_mesh%nodes%comm%gather(loc_coords,coords)

    if( myproc .eq. 0) then
      open(11,file=filename,access='sequential',status='unknown')
      open(21,file=rtable,access='sequential',status='unknown')

      read(21,*)
      read(21,*)ndgl
      allocate(nloen(ndgl))
      do jlat=1,ndgl
        read(21,*) idum,nloen(jlat)
      enddo
      nb_nodes = sum(nloen)+ndgl ! plus periodic points
      allocate( plot_idx(nb_nodes) )
      inode = 1
      jglb = 1
      do jlat=1,ndgl
        inode_grwch = inode
        do jnode=1,nloen(jlat)
          plot_idx(inode) = jglb
          inode = inode+1
          jglb = jglb+1
        end do
        plot_idx(inode) = plot_idx(inode_grwch)
        inode = inode+1
      enddo
      write(11,*)g%fields%time
      write(11,*)0
      write(11,'(A)')'point X Y p'


      do jnode=1,nb_nodes
        inode = plot_idx(jnode)
        write(11,*) jnode, coords(inode,XX), coords(inode,YY), &
          & D(inode), Q(inode,XX)/D(inode), Q(inode,YY)/D(inode), D(inode)
      enddo
      close(21)
    end if
  end subroutine write_results_joanna

end module joanna_module
