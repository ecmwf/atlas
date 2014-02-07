
! module for testing creation of a mesh
module gmsh_module
  use parallel_module, only: myproc, nproc, parallel_barrier, Comm_type
  use datastruct_module
  use datastruct
  use common_module, only : jprw,log_info,XX,YY,ZZ, str
  implicit none

  real(kind=jprw), parameter :: pi = acos(-1.)

contains

  subroutine read_gmsh(filename,dstruct)
    character(len=*), intent(in) :: filename
    type(DataStructure_type), intent(inout) :: dstruct
    type(Field_type) :: field
    integer , pointer :: bounds(:)
    integer :: jnode, jedge, ip1, ip2, ibound
    integer , allocatable :: tmp(:)
    real(kind=jprw), pointer :: edge_centroids(:,:)
    dstruct%mesh = datastruct_read_gmsh(filename)
    call datastruct_build_periodic_boundaries(dstruct%mesh);
    call datastruct_build_edges(dstruct%mesh);
    call datastruct_build_dual_mesh(dstruct%mesh);

    dstruct%functionspace_nodes_2d = dstruct%mesh%function_space("nodes_2d")
    dstruct%functionspace_edges_2d = dstruct%mesh%function_space("edges")

    dstruct%time = 0.
    dstruct%time_step = 0
    dstruct%fields = new_FieldSet("fields")
    dstruct%output_fields = new_FieldSet("output_fields")

    bounds => dstruct%functionspace_nodes_2d%bounds()
    dstruct%nb_nodes = bounds(2)

    bounds => dstruct%functionspace_edges_2d%bounds()
    dstruct%nb_edges = bounds(2)

    field = dstruct%functionspace_edges_2d%field("nodes")
    call field%access_data( dstruct%edges )

    field = dstruct%functionspace_nodes_2d%field("proc")
    call field%access_data( dstruct%nodes_proc )

    field = dstruct%functionspace_nodes_2d%field("glb_idx")
    call field%access_data( dstruct%nodes_glb_idx )

    field = dstruct%functionspace_edges_2d%field("proc")
    call field%access_data( dstruct%edges_proc )

    field = dstruct%functionspace_edges_2d%field("glb_idx")
    call field%access_data( dstruct%edges_glb_idx )

    call dstruct%fields%add_field( dstruct%functionspace_nodes_2d%field("coordinates") )
    call dstruct%fields%add_field( dstruct%functionspace_nodes_2d%field("dual_volumes") )
    call dstruct%fields%add_field( dstruct%functionspace_edges_2d%field("dual_normals") )


    call dstruct%nodes_comm%setup( dstruct%nodes_proc, dstruct%nodes_glb_idx )


    ! Connectivity setup in C++ needs to be changed to fortran indexing starting from 1
    do jedge=1,dstruct%nb_edges
      dstruct%edges(:,jedge) = dstruct%edges(:,jedge) + 1
    enddo

    dstruct%glb_nb_edges = dstruct%nb_edges

    allocate(dstruct%nb_neighbours(dstruct%nb_nodes))
    dstruct%nb_neighbours(:) = 0
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(1,jedge)
      ip2 = dstruct%edges(2,jedge)
      dstruct%nb_neighbours(ip1) = dstruct%nb_neighbours(ip1)+1
      dstruct%nb_neighbours(ip2) = dstruct%nb_neighbours(ip2)+1
    enddo
    dstruct%max_nb_neighbours = maxval(dstruct%nb_neighbours(:))
    ibound = 0
    do jnode = 1,dstruct%nb_nodes
      if( dstruct%nb_neighbours(jnode) == 1) ibound=ibound+1
    enddo

    allocate(dstruct%neighbours(dstruct%max_nb_neighbours,dstruct%nb_nodes))
    allocate(dstruct%my_edges(dstruct%max_nb_neighbours,dstruct%nb_nodes))
    allocate(dstruct%sign(dstruct%max_nb_neighbours,dstruct%nb_nodes))
    dstruct%nb_neighbours(:) = 0
    allocate( tmp(dstruct%nb_edges) )

    field = dstruct%functionspace_edges_2d%field("centroids")
    call field%access_data( edge_centroids )

    dstruct%nb_pole_edges = 0
    do jedge = 1,dstruct%nb_edges
      ip1 = dstruct%edges(1,jedge)
      ip2 = dstruct%edges(2,jedge)
      dstruct%nb_neighbours(ip1) = dstruct%nb_neighbours(ip1)+1
      dstruct%neighbours(dstruct%nb_neighbours(ip1),ip1) = ip2
      dstruct%my_edges(dstruct%nb_neighbours(ip1),ip1) = jedge
      dstruct%sign(dstruct%nb_neighbours(ip1),ip1) = 1
      dstruct%nb_neighbours(ip2) = dstruct%nb_neighbours(ip2)+1
      dstruct%neighbours(dstruct%nb_neighbours(ip2),ip2) = ip1
      dstruct%my_edges(dstruct%nb_neighbours(ip2),ip2) = jedge
      dstruct%sign(dstruct%nb_neighbours(ip2),ip2) = -1
      if( abs ( abs(edge_centroids(YY,jedge)) - 0.5_jprw*pi ) < 1.e-6_jprw ) then
        dstruct%nb_pole_edges = dstruct%nb_pole_edges + 1
        tmp(dstruct%nb_pole_edges) = jedge
      endif
    enddo
    allocate(dstruct%pole_edges(dstruct%nb_pole_edges))
    dstruct%pole_edges(:) = tmp(1:dstruct%nb_pole_edges)
    write(0,*) "nb_pole_edges = " , dstruct%nb_pole_edges


    ! Pole edges still missing
    !In C++ compute centroids of pole edges as y=pi/2 (or something)
    !Identify pole edges here, and allocate it




  end subroutine read_gmsh

  subroutine write_gmsh_nodal_field_2d(field, comm, dstruct)
    type(Field_type), intent(in) :: field
    type(Comm_type), intent(in) :: comm
    class(DataStructure_type), intent(in) :: dstruct
    
    real(kind=jprw), allocatable :: glb_field(:,:)
    integer :: jnode, glb_rows

    glb_rows = comm%glb_field_size()

    if ( nproc > 1 ) then
      call comm%gather( field%data2(), glb_field )
    else
      allocate(glb_field(field%nb_vars(),dstruct%nb_nodes))
      glb_field = field%data2()
    endif

    if (myproc .eq. 0) then
      write(51,'(A)')"$NodeData"
      write(51,*) 1                     ! one string tag:
      write(51,'(A)') '"'//field%name()//'"'      ! the name of the view ("A scalar view")
      write(51,*) 1                     ! one real(kind=jprw) tag:
      write(51,*) dstruct%time          ! the time value (0.0)
      write(51,*) 3                     ! three integer tags:
      write(51,*) dstruct%time_step     ! the time step (0; time steps always start at 0)
      if (field%nb_vars() == 1) write(51,*) 1    ! 1-component (scalar) field
      if (field%nb_vars() == 2) write(51,*) 3    ! 3-component (vector) field
      write(51,*) size(glb_field,2)          ! number of associated nodal values

      if (field%nb_vars() == 1) then
        do jnode=1,size(glb_field,2)
          write(51,'(1I8,E18.10)') jnode, glb_field(1,jnode)
        enddo
      endif
      if (field%nb_vars() == 2) then
        do jnode=1,glb_rows
          write(51,'(1I8,E18.10,E18.10,F3.0)') jnode, glb_field(1,jnode), glb_field(2,jnode), 0.0
        enddo
      endif
      write(51,'(A)')"$EndNodeData"
    end if
  end subroutine write_gmsh_nodal_field_2d
  
  subroutine write_gmsh_mesh_2d(dstruct,filename)
    class(DataStructure_type), intent(in) :: dstruct
    character(len=*), intent(in) :: filename
    character(len=1024) :: procfile
    integer :: jedge, jnode
    real(kind=jprw) :: r, phi, theta, x, y, z, pi
    real(kind=jprw), pointer :: coords(:,:)
    pi = acos(-1._jprw)

    coords => vector_field_2d("coordinates",dstruct)

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
      phi   = coords(XX,jnode)
      theta = -(coords(YY,jnode)+pi/2.)
      x = r*sin(theta)*cos(phi)
      y = r*sin(theta)*sin(phi)
      z = r*cos(theta)

      write(50,'(1I8,F18.10,F18.10,F18.10)')  dstruct%nodes_glb_idx(jnode), x,y,z
      !write(50,'(1I8,F18.10,F18.10,F18.10)') dstruct%nodes_glb_idx(jnode), &
      !  & coords(XX,jnode), coords(YY,jnode), 0.
    enddo
    write(50,'(A)')"$EndNodes"
    write(50,'(A)')"$Elements"
    write(50,*) dstruct%nb_edges
    do jedge=1, dstruct%nb_edges
      ! element-number  type(1=lineP1)  nb_tags(=3)  tag1(=physical-group)  tag2(=elementary-group) tag4(=nb_partitions) tag3(=partition) [nodes]
      write(50,'(1I8,5I2,I3, 2I8)')  dstruct%edges_glb_idx(jedge), 1, 4, 1, 1, 1, dstruct%edges_proc(jedge), &
        & dstruct%nodes_glb_idx( dstruct%edges(1,jedge) ), dstruct%nodes_glb_idx( dstruct%edges(2,jedge) )
    enddo
    write(50,'(A)')"$EndElements"
    close(50)
  end subroutine write_gmsh_mesh_2d

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
      if (function_space%name() == "nodes_2d") then
        call log_info("writing 2d field "//trim(field%name()))
        call write_gmsh_nodal_field_2d(field,dstruct%nodes_comm,dstruct)
      else if (function_space%name() == "nodes_3d") then
        call log_info("writing 3d field "//trim(field%name()))
        call write_gmsh_nodal_field_3d(field,dstruct%nodes_comm,dstruct)
      end if
    end do
      
    if (myproc .eq. 0) then
      close(51)
    end if
  end subroutine write_gmsh_fields


  subroutine write_gmsh_nodal_field_3d(field, comm, dstruct)
    type(Field_type), intent(in) :: field
    type(Comm_type), intent(in) :: comm
    type(DataStructure_type), intent(in) :: dstruct
    
    real(kind=jprw), pointer :: loc_field(:,:,:)
    real(kind=jprw), allocatable :: glb_field(:,:,:)
    integer :: jnode, jlev, nb_levels

    loc_field => field%data3()

    call comm%gather( field%data3(), glb_field )
    nb_levels = dstruct%nb_levels

  do jlev=1,nb_levels

    if (myproc .eq. 0) then
      write(51,'(A)')"$NodeData"
      write(51,*) 1                     ! one string tag:
      write(51,'(A)') '"'//field%name()//trim(str(jlev,'(I3)'))//'"'      ! the name of the view ("A scalar view")
      write(51,*) 1                     ! one real(kind=jprw) tag:
      write(51,*) dstruct%time          ! the time value (0.0)
      write(51,*) 3                     ! three integer tags:
      write(51,*) dstruct%time_step     ! the time step (0; time steps always start at 0)
      if (field%nb_vars() == 1) write(51,*) 1    ! 1-component (scalar) field
      if (field%nb_vars() == 2) write(51,*) 3    ! 3-component (vector) field
      if (field%nb_vars() == 3) write(51,*) 3    ! 3-component (vector) field
      write(51,*) size(glb_field,3)          ! number of associated nodal values

      !do jlev=1,nb_levels
        do jnode=1, dstruct%glb_nb_nodes
          if (field%nb_vars() == 1) then
            write(51,'(1I8,E18.10)') jnode, glb_field(1,jlev,jnode)
          endif
          if (field%nb_vars() == 2) then
            write(51,'(1I8,E18.10,E18.10,F3.0)') jnode, glb_field(1,jlev,jnode), glb_field(2,jlev,jnode), 0.0
          endif
          if (field%nb_vars() == 3) then
            write(51,'(1I8,E18.10,E18.10,F3.0)') jnode, glb_field(1,jlev,jnode), glb_field(2,jlev,jnode), glb_field(3,jlev,jnode)
          endif
        end do
      !end do
      write(51,'(A)')"$EndNodeData"
    end if
  end do
  end subroutine write_gmsh_nodal_field_3d



  subroutine write_gmsh_mesh_3d(dstruct,filename)
    class(DataStructure_type), intent(in) :: dstruct
    character(len=*), intent(in) :: filename
    character(len=1024) :: procfile
    integer :: jedge, jnode, jlev, gidn, gide, gidp1, gidp2
    real(kind=jprw) :: r, phi, theta, x, y, z, pi
    real(kind=jprw), pointer :: coords(:,:)
    pi = acos(-1._jprw)

    coords => vector_field_2d("coordinates",dstruct)

    call log_info( "Writing Gmsh file "//trim(filename) )
    
    call parallel_barrier()

    write(procfile,'(A,A,I2.2)') trim(filename),'.P',myproc
    open(50,file=trim(procfile),access='sequential',status='REPLACE')
    
    write(50,'(A)')"$MeshFormat"
    write(50,'(A)')"2.2 0 8"
    write(50,'(A)')"$EndMeshFormat"
    write(50,'(A)')"$Nodes"
    write(50,*) dstruct%nb_nodes * dstruct%nb_levels

    do jlev=1, dstruct%nb_levels
      r     = 6371.! + 6371.*.1*(jlev-1._jprw)
      do jnode=1, dstruct%nb_nodes
        gidn = (jlev-1)*dstruct%glb_nb_nodes + dstruct%nodes_glb_idx(jnode)
        phi   = coords(XX,jnode)
        theta = -(coords(YY,jnode)+pi/2.)
        x = r*sin(theta)*cos(phi)
        y = r*sin(theta)*sin(phi)
        z = r*cos(theta)
        write(50,'(1I8,F18.10,F18.10,F18.10)') gidn, x,y,z
        !write(50,'(1I8,F18.10,F18.10,F18.10)')  mesh%nodes%glb_idx(jnode), &
        !  & mesh%nodes_coordinates(jnode,1),mesh%nodes_coordinates(jnode,2), 0.
      end do
    enddo
    write(50,'(A)')"$EndNodes"
    write(50,'(A)')"$Elements"
    write(50,*) dstruct%nb_edges * dstruct%nb_levels
    do jlev=1, dstruct%nb_levels
      do jedge=1, dstruct%nb_edges
        gide  = (jlev-1)*dstruct%glb_nb_edges + dstruct%edges_glb_idx(jedge)
        gidp1 = (jlev-1)*dstruct%glb_nb_nodes + dstruct%nodes_glb_idx(dstruct%edges(jedge,1))
        gidp2 = (jlev-1)*dstruct%glb_nb_nodes + dstruct%nodes_glb_idx(dstruct%edges(jedge,2))

        ! element-number  type(1=lineP1)  nb_tags(=4)  tag1(=physical-group)  tag2(=elementary-group) tag3(=nb_partitions) tag4(=partition) [nodes]
        write(50,'(1I8,5I2,I3, 2I8)') gide, 1, 4, jlev, jlev, 1, dstruct%edges_proc(jedge), &
          & gidp1, gidp2
      enddo
    end do
    write(50,'(A)')"$EndElements"
    close(50)
  end subroutine write_gmsh_mesh_3d

end module gmsh_module
