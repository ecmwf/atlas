
! module for testing creation of a Grid
module example_grids
  use grid_module
  implicit none
contains
    
  ! read_quads(grid)
  ! -----------------
  ! Create a example grid containing 3 quadrilaterals
  !
  !  4------3------6------8
  !  |      |      |      |
  !  |  E1  |  E2  |  E3  |
  !  |      |      |      |
  !  1------2------5------7
  subroutine read_quads(g)
    class(Grid_class), intent(inout) :: g
      
    ! Specify element type
    call g%init("LagrangeP1_Quad2D")
    
    ! Specify number of cells and nodes
    g%nb_elems =  3
    g%nb_nodes =  8
    g%nb_faces = 10
    
    ! Fill in nodes and cells
    allocate(g%nodes(g%nb_nodes,g%dimension))
    g%nodes(1,:) = [0,0] 
    g%nodes(2,:) = [1,0] 
    g%nodes(3,:) = [1,1] 
    g%nodes(4,:) = [0,1]
    g%nodes(5,:) = [2,0] 
    g%nodes(6,:) = [2,1] 
    g%nodes(7,:) = [3,0] 
    g%nodes(8,:) = [3,1]

    allocate(g%cells(g%nb_elems,g%cell%nb_nodes))
    g%cells(1,:) = [1,2,3,4]
    g%cells(2,:) = [2,5,6,3]
    g%cells(3,:) = [5,7,8,6]

    allocate(g%faces(g%nb_faces,2)) ! this 2 should be nb_nodes_per_face
    g%faces( 1,:) = [1,2]
    g%faces( 2,:) = [2,3]
    g%faces( 3,:) = [3,4]
    g%faces( 4,:) = [4,1]
    g%faces( 5,:) = [2,5]
    g%faces( 6,:) = [5,6]
    g%faces( 7,:) = [6,3]
    g%faces( 8,:) = [5,7]
    g%faces( 9,:) = [7,8]
    g%faces(10,:) = [8,6]

  end subroutine read_quads
  
  ! read_quads(grid)
  ! -----------------
  ! Create a example grid containing 4 triangles
  !
  !  3------4------6
  !  | \  E2| \  E4|
  !  |   \  |   \  |
  !  |E1   \|E3   \|
  !  1------2------5
  subroutine read_triags(g)
    class(Grid_class), intent(inout) :: g
      
    ! Specify element type
    call g%init("LagrangeP1_Triag2D")
    
    ! Specify number of cells and nodes
    g%nb_elems = 4
    g%nb_nodes = 6
    g%nb_faces = 9
    
    ! Fill in nodes and cells
    allocate(g%nodes(g%nb_nodes,g%dimension))
    g%nodes(1,:) = [0,0] 
    g%nodes(2,:) = [2,0] 
    g%nodes(3,:) = [0,2]
    g%nodes(4,:) = [2,2] 
    g%nodes(5,:) = [4,0] 
    g%nodes(6,:) = [4,2] 

    allocate(g%cells(g%nb_elems,g%cell%nb_nodes))
    g%cells(1,:) = [1,2,3]
    g%cells(2,:) = [2,4,3]
    g%cells(3,:) = [2,5,4]
    g%cells(4,:) = [5,6,4]

    allocate(g%faces(g%nb_faces,2)) ! this 2 should be nb_nodes_per_face
    g%faces( 1,:) = [1,2]
    g%faces( 2,:) = [2,3]
    g%faces( 3,:) = [3,1]
    g%faces( 4,:) = [2,3]
    g%faces( 5,:) = [4,3]
    g%faces( 6,:) = [2,5]
    g%faces( 7,:) = [5,4]
    g%faces( 8,:) = [5,6]
    g%faces( 9,:) = [6,4]

  end subroutine read_triags
  
end module example_grids

! =============================================================================
! =============================================================================
! =============================================================================

! Main program
program main
  use grid_module
  use example_grids
  implicit none
  
  ! Declarations
  ! ------------
  class(Grid_class)                        :: g
  class(FunctionSpace),  pointer    :: cellcentred  
  class(FunctionSpace),  pointer    :: vertices
  class(FunctionSpace),  pointer    :: faces

  type(Field),           pointer    :: X, DxDx
  type(Field),           pointer    :: T
  type(Field),           pointer    :: State
  type(Field),           pointer    :: S

  integer                           :: n,e
  real, dimension(:),   allocatable :: local_coord
  real, dimension(:),   allocatable :: physical_coord

  real, dimension(:,:), allocatable :: cell_coords
  real, dimension(:),   allocatable :: interpolate
  real, dimension(:,:), allocatable :: local_gradient
  real, dimension(:,:), allocatable :: physical_gradient
  real, dimension(:,:), allocatable :: Jinv
  real, dimension(:,:), allocatable :: grad

 
  ! Execution
  ! ---------
  
  ! Read grid from example_grids module
 
  !call read_quads(g)
  call read_triags(g)

  faces => new_FaceFunctionSpace("faces", "LagrangeP0", g)

  ! Setup a discontinuous function space of P0 cells ( CELLCENTRED )
  cellcentred => new_DiscontinuousFunctionSpace("cellcentred","LagrangeP0", g)
  
  ! Setup a continuous function space of P1 cells ( VERTICES )
  vertices => new_ContinuousFunctionSpace("vertices","LagrangeP1", g)

  ! Create a scalar field for temperature in the function spaces
  X     => cellcentred%add_vector_field("X")
  DxDx  => cellcentred%add_scalar_field("DxDx")
  T     => vertices%add_scalar_field("Temperature")
  State => vertices%add_array_field("State",4)
  S     => faces%add_vector_field("face_normal")

  ! Allocations
  allocate( physical_coord(g%cell%dimension) )
  allocate( local_coord(g%cell%dimensionality) )
  allocate( cell_coords(g%cell%nb_nodes, g%cell%dimension) )
  allocate( Jinv(g%cell%dimension,g%cell%dimensionality) )
  allocate( interpolate(vertices%sf%nb_nodes) )
  allocate( local_gradient(vertices%sf%dimensionality,vertices%sf%nb_nodes) )
  allocate( physical_gradient(g%dimension,vertices%sf%nb_nodes) )
  allocate( grad(g%dimension,g%dimension) )

  ! Use iso-parametric information to interpolate and compute gradient in element local coords
  call g%cell%sf%values( cellcentred%sf%local_coords(1,:), interpolate )
  call g%cell%sf%grad_values( cellcentred%sf%local_coords(1,:), local_gradient )

  write(0,*) "interpolate ",interpolate
  write(0,*) "local_gradient ",local_gradient

  ! Loop over cells
  do e=1,g%nb_elems

    call g%cell_coords(e, cell_coords)
    physical_coord = matmul(interpolate, cell_coords)
    X%array(e,:) = physical_coord

    call g%cell%jacobian_inverse( cellcentred%sf%local_coords(1,:), cell_coords, Jinv )
    physical_gradient = matmul(Jinv,local_gradient)
    grad = matmul(physical_gradient, cell_coords) ! [ dxdx dydx; dxdy dydy ]
    DxDx%array(e,:) = grad(1,1)

  end do

  ! Do some science with the temperature field
  do n=1,T%size
     T%array(n,1) = 20+n
  end do

  ! Print
  ! -----
  write(0,*) X%name," (",X%size,",",X%cols,") = "
  do n=1,X%size
    write(0,*) X%array(n,:)
  end do
  write(0,*) DxDx%name," (",DxDx%size,",",DxDx%cols,") = "
  do n=1,DxDx%size
    write(0,*) DxDx%array(n,:)
  end do
  write(0,*) T%name," (",T%size,",",T%cols,") = "
  do n=1,T%size
    write(0,*) T%array(n,:)
  end do
  write(0,*) State%name," (",State%size,",",State%cols,") = "
  do n=1,State%size
    write(0,*) State%array(n,:)
  end do
    write(0,*) S%name," (",S%size,",",S%cols,") = "
  do n=1,S%size
    write(0,*) S%array(n,:)
  end do
  
  ! Destruction
  ! -----------
  ! Recursively deallocate the grid, functionspaces, fields, ...
  call g%destruct()
  deallocate(physical_coord)
  deallocate(local_coord)
  deallocate(cell_coords)
  deallocate(Jinv)
  deallocate(interpolate)
  deallocate(local_gradient)
  deallocate(physical_gradient)
  deallocate(grad)
end program main