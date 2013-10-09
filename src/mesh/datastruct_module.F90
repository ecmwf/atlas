
module datastruct_module

  use grid_module

  implicit none
  private
  public :: create_mesh
  public :: create_scalar_field_in_nodes
  public :: create_vector_field_in_nodes
  public :: create_scalar_field_in_edges
  public :: create_vector_field_in_edges
  public :: scalar_field, vector_field
  public :: geom

  type, public, extends(Grid_class) :: Mesh_class
  end type Mesh_class

  type, public :: DataStructure_type
  private
  ! PRIVATE part
  ! ------------
    type(Mesh_class), public     :: internal_mesh
    type(State_class), public    :: fields
    class(FunctionSpace_class), pointer     :: functionspace_nodes
    class(FunctionSpace_class), pointer     :: functionspace_edges
    type(Field_class)                       :: field_coordinates
    type(Field_class)                       :: field_dual_volumes
    type(Field_class)                       :: field_dual_normals

  ! PUBLIC part
  ! -----------
    integer, public                                  :: nb_nodes=0
    integer, public                                  :: nb_edges=0
    integer, public                                  :: nb_internal_edges=0
    integer, public                                  :: nb_pole_edges=0
    integer, public                                  :: nb_ghost_nodes=0

    integer,         dimension(:,:), allocatable, public :: edges
    integer,         dimension(:),   allocatable, public :: internal_edges
    integer,         dimension(:),   allocatable, public :: pole_edges
    integer,         dimension(:,:), allocatable, public :: ghost_nodes
    real(kind=jprb), dimension(:,:), allocatable, public :: coordinates  ! Defined in nodes
    real(kind=jprb), dimension(:),   allocatable, public :: dual_volumes ! Defined in nodes
    real(kind=jprb), dimension(:,:), allocatable, public :: dual_normals ! Defined in edges

  end type DataStructure_type

  type(DataStructure_type) :: geom

contains
  
  subroutine create_mesh(nb_nodes, nb_edges, geom_)
    implicit none
    integer, intent(in) :: nb_nodes
    integer, intent(in) :: nb_edges
    type(DataStructure_type), intent(inout), target :: geom_

    call geom_%internal_mesh%init("LagrangeP1_Triag2D")
    allocate(geom_%internal_mesh%nodes_coordinates(nb_nodes,2))
    call geom_%fields%init("geometry_fields")

    geom_%nb_nodes = nb_nodes
    geom_%internal_mesh%nb_nodes = nb_nodes
    allocate( ContinuousFunctionSpace :: geom_%functionspace_nodes )
    call geom_%functionspace_nodes%init("nodes", "LagrangeP1", geom_%internal_mesh)
    
    !call geom_%field_coordinates%init("coordinates",geom_%functionspace_nodes,2)
    !geom_%coordinates => geom_%field_coordinates%array
    allocate( geom_%coordinates(nb_nodes,2) )
    !call geom_%field_dual_volumes%init("dual_volumes",geom_%functionspace_nodes,1)
    !geom_%dual_volumes => geom_%field_dual_volumes%array(:,1)
    allocate( geom_%dual_volumes(nb_nodes) )
    
    call geom_%internal_mesh%add_function_space(geom_%functionspace_nodes)

    geom_%nb_edges = nb_edges
    geom_%internal_mesh%nb_faces = nb_edges
    allocate( FaceFunctionSpace :: geom_%functionspace_edges )
    call geom_%functionspace_edges%init("edges", "LagrangeP0", geom_%internal_mesh)
    !call geom_%field_dual_normals%init("dual_normals", geom_%functionspace_edges,2)
    !geom_%dual_normals => geom_%field_dual_normals%array
    allocate( geom_%dual_normals(nb_edges,2) )
    call geom_%internal_mesh%add_function_space(geom_%functionspace_edges)
    allocate( geom_%internal_mesh%faces(geom_%nb_edges,2) )
    allocate( geom_%edges(geom_%nb_edges,2) )
    !geom_%edges => geom_%internal_mesh%faces
  end subroutine create_mesh

  subroutine create_scalar_field_in_nodes(name, geom_)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom_
    call geom_%fields%add_field(name, geom_%functionspace_nodes,1)
  end subroutine create_scalar_field_in_nodes

  subroutine create_vector_field_in_nodes(name, geom_)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom_
    call geom_%fields%add_field(name, geom_%functionspace_nodes,2)
  end subroutine create_vector_field_in_nodes

  subroutine create_scalar_field_in_edges(name, geom_)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom_
    call geom_%fields%add_field(name, geom_%functionspace_edges,1)
  end subroutine create_scalar_field_in_edges

  subroutine create_vector_field_in_edges(name, geom_)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom_
    call geom_%fields%add_field(name, geom_%functionspace_edges,2)
  end subroutine create_vector_field_in_edges

  function scalar_field(name, geom_) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom_
    type(Field_class), pointer :: field
    real(kind=jprb), dimension(:), pointer :: array
    field => geom_%fields%field(name)
    array => field%array(:,1)
  end function scalar_field

  function vector_field(name, geom_) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom_
    type(Field_class), pointer :: field
    real(kind=jprb), dimension(:,:), pointer :: array
    field => geom_%fields%field(name)
    array => field%array
  end function vector_field

end module datastruct_module