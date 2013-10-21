module datastruct_module

  use common_module, only: jprb
  use grid_module, only : Grid_class, State_class, Field_class, FaceFunctionSpace, ContinuousFunctionSpace, FunctionSpace_class
  use parallel_module

  implicit none
  private
  public :: create_mesh
  public :: create_scalar_field_in_nodes
  public :: create_vector_field_in_nodes
  public :: create_scalar_field_in_edges
  public :: create_vector_field_in_edges
  public :: scalar_field, vector_field
  public :: synchronise

  interface  synchronise
    module procedure synchronise_field
    module procedure synchronise_array_rank1
    module procedure synchronise_array_rank2
  end interface synchronise
!  type, public, extends(Grid_class) :: Mesh_class
!  end type Mesh_class

  type, public :: DataStructure_type
  private
  ! PRIVATE part
  ! ------------
    type(Grid_class), public     :: internal_mesh
    type(State_class), public    :: fields
    class(FunctionSpace_class), pointer     :: functionspace_nodes
    class(FunctionSpace_class), pointer     :: functionspace_edges

  ! PUBLIC part
  ! -----------
    integer, public                                  :: glb_nb_nodes=0
    integer, public                                  :: glb_nb_edges=0
    integer, public                                  :: nb_nodes=0
    integer, public                                  :: nb_edges=0
    integer, public                                  :: nb_pole_edges=0
    integer, public                                  :: nb_ghost_nodes=0

    integer, dimension(:,:), pointer, public     :: edges
    integer, dimension(:),   allocatable, public :: pole_edges
    integer, dimension(:,:), allocatable, public :: ghost_nodes

  end type DataStructure_type


contains
  
  subroutine create_mesh(nb_nodes, nb_edges, proc, glb_idx, geom)
    implicit none
    integer, intent(in) :: nb_nodes
    integer, intent(in) :: nb_edges
    integer, intent(in) :: proc(:)
    integer, intent(in) :: glb_idx(:)
    type(DataStructure_type), intent(inout), target :: geom

    call geom%internal_mesh%init("LagrangeP1_Triag2D")
    allocate(geom%internal_mesh%nodes_coordinates(nb_nodes,2))
    call geom%fields%init("geometry_fields")

    geom%nb_nodes = nb_nodes
    geom%nb_edges = nb_edges
    geom%internal_mesh%nb_nodes = nb_nodes
    geom%internal_mesh%nb_faces = nb_edges
    allocate( geom%internal_mesh%faces(geom%nb_edges,2) )
    geom%edges => geom%internal_mesh%faces

    ! Create 1 function space for nodes, (DONT WORRY ABOUT DETAILS)
    allocate( ContinuousFunctionSpace :: geom%functionspace_nodes )
    geom%internal_mesh%nodes =>  geom%functionspace_nodes 
    call geom%functionspace_nodes%init("nodes", "LagrangeP1", geom%internal_mesh)
    call geom%internal_mesh%add_function_space(geom%functionspace_nodes)

    ! Create 1 function space for edges, (DONT WORRY ABOUT DETAILS)
    allocate( FaceFunctionSpace :: geom%functionspace_edges )
    call geom%functionspace_edges%init("edges", "LagrangeP0", geom%internal_mesh)
    call geom%internal_mesh%add_function_space(geom%functionspace_edges)

    call create_vector_field_in_nodes("coordinates",geom)
    call create_scalar_field_in_nodes("dual_volumes",geom)
    call create_vector_field_in_edges("dual_normals", geom)

    allocate( geom%functionspace_nodes%proc( size(proc) ) )
    allocate( geom%functionspace_nodes%glb_idx( size(glb_idx ) ) )
    geom%functionspace_nodes%proc(:) = proc(:)
    geom%functionspace_nodes%glb_idx(:) = glb_idx(:)
    call geom%functionspace_nodes%comm%setup( geom%functionspace_nodes%proc, geom%functionspace_nodes%glb_idx )
  end subroutine create_mesh

  subroutine create_scalar_field_in_nodes(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%fields%add_field(name, geom%functionspace_nodes,1)
  end subroutine create_scalar_field_in_nodes

  subroutine create_vector_field_in_nodes(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%fields%add_field(name, geom%functionspace_nodes,2)
  end subroutine create_vector_field_in_nodes

  subroutine create_scalar_field_in_edges(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%fields%add_field(name, geom%functionspace_edges,1)
  end subroutine create_scalar_field_in_edges

  subroutine create_vector_field_in_edges(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%fields%add_field(name, geom%functionspace_edges,2)
  end subroutine create_vector_field_in_edges

  function scalar_field(name, geom) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    type(Field_class), pointer :: field
    real(kind=jprb), dimension(:), pointer :: array
    field => geom%fields%field(name)
    array => field%array(:,1)
  end function scalar_field

  function vector_field(name, geom) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    type(Field_class), pointer :: field
    real(kind=jprb), dimension(:,:), pointer :: array
    field => geom%fields%field(name)
    array => field%array
  end function vector_field

  subroutine synchronise_field(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    type(Field_class), pointer :: field
    field => geom%fields%field(name)
    call field%function_space%comm%synchronise( field%array )
  end subroutine synchronise_field

  subroutine synchronise_array_rank1(array, geom)
    implicit none
    real(kind=jprb), dimension(:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: geom
    call geom%functionspace_nodes%comm%synchronise( array )
  end subroutine synchronise_array_rank1

  subroutine synchronise_array_rank2(array, geom)
    implicit none
    real(kind=jprb), dimension(:,:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: geom
    call geom%functionspace_nodes%comm%synchronise( array )
  end subroutine synchronise_array_rank2

end module datastruct_module