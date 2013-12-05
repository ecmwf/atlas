module datastruct_module

  use common_module, only: jprw
  use datastruct, only : Mesh_type, FieldSet_type, Field_type, &
      & new_FunctionSpace, FunctionSpace_type, new_Mesh, new_FieldSet
  use parallel_module

  implicit none
  private
  public :: create_mesh
  public :: create_scalar_field_in_nodes
  public :: create_vector_field_in_nodes
  public :: create_scalar_field_in_edges
  public :: create_vector_field_in_edges
  public :: scalar_field, vector_field
  public :: mark_output
  public :: synchronise, gather

  interface  synchronise
    module procedure synchronise_array_rank1
    module procedure synchronise_array_rank2
  end interface synchronise

  interface  gather
    module procedure gather_array_rank1
    module procedure gather_array_rank2
  end interface gather
!  type, public, extends(Mesh_type) :: Mesh_class
!  end type Mesh_class

  type, public :: DataStructure_type
  private
  ! PRIVATE part
  ! ------------
    type(Mesh_type)     :: mesh
    type(FieldSet_type) :: fields
    type(FunctionSpace_type) :: functionspace_nodes
    type(FunctionSpace_type) :: functionspace_edges
    type(Comm_type), public :: nodes_comm

    type(FieldSet_type), public :: output_fields

    real(kind=jprw), public :: time

  ! PUBLIC part
  ! -----------
    integer,                       public :: glb_nb_nodes=0
    integer,                       public :: glb_nb_edges=0
    integer,                       public :: nb_nodes=0
    integer,                       public :: nb_edges=0
    integer,                       public :: nb_pole_edges=0
    integer,                       public :: nb_ghost_nodes=0

    integer,          pointer,     public :: edges(:,:)
    integer,          allocatable, public :: pole_edges(:)
    integer,          allocatable, public :: ghost_nodes(:,:)
    integer,          allocatable, public :: nb_neighbours(:)
    integer,          allocatable, public :: neighbours(:,:)
    integer,          allocatable, public :: my_edges(:,:)
    real(kind=jprw),  allocatable, public :: sign(:,:)
    integer,                       public :: max_nb_neighbours        

    integer,          allocatable, public :: nodes_proc(:)
    integer,          allocatable, public :: nodes_glb_idx(:)
    integer,          allocatable, public :: edges_proc(:)
    integer,          allocatable, public :: edges_glb_idx(:)

  end type DataStructure_type


contains
  
  subroutine create_mesh(nb_nodes, nb_edges, proc, glb_idx, geom)
    implicit none
    integer, intent(in) :: nb_nodes
    integer, intent(in) :: nb_edges
    integer, intent(in) :: proc(:)
    integer, intent(in) :: glb_idx(:)
    type(DataStructure_type), intent(inout), target :: geom

    geom%mesh = new_Mesh()
    geom%fields = new_FieldSet("fields")
    geom%output_fields = new_FieldSet("output_fields")

    call geom%mesh%add_function_space( new_FunctionSpace("nodes", "P1-C", nb_nodes) )
    geom%functionspace_nodes = geom%mesh%function_space("nodes")
    geom%nb_nodes = nb_nodes

    call geom%mesh%add_function_space( new_FunctionSpace("edges", "P0-D", nb_edges) )
    geom%functionspace_edges = geom%mesh%function_space("edges")
    geom%nb_edges = nb_edges

    allocate( geom%edges(geom%nb_edges,2) )

    call create_vector_field_in_nodes("coordinates",geom)
    call create_scalar_field_in_nodes("dual_volumes",geom)
    call create_vector_field_in_edges("dual_normals", geom)

    allocate( geom%nodes_proc( nb_nodes ) ) 
    allocate( geom%nodes_glb_idx( nb_nodes ) ) 
    geom%nodes_proc(:) = proc(:)
    geom%nodes_glb_idx(:) = glb_idx(:)
    call geom%nodes_comm%setup( proc, glb_idx )

    allocate( geom%edges_proc( nb_edges ) ) 
    allocate( geom%edges_glb_idx( nb_edges ) ) 

  end subroutine create_mesh


  subroutine create_scalar_field_in_nodes(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%functionspace_nodes%create_field(name,1)
    call geom%fields%add_field( geom%functionspace_nodes%field(name) )
  end subroutine create_scalar_field_in_nodes

  subroutine create_vector_field_in_nodes(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%functionspace_nodes%create_field(name,2)
    call geom%fields%add_field( geom%functionspace_nodes%field(name) )
  end subroutine create_vector_field_in_nodes

  subroutine create_scalar_field_in_edges(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%functionspace_edges%create_field(name,1)
    call geom%fields%add_field( geom%functionspace_edges%field(name) )
  end subroutine create_scalar_field_in_edges

  subroutine create_vector_field_in_edges(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(inout) :: geom
    call geom%functionspace_edges%create_field(name,2)
    call geom%fields%add_field( geom%functionspace_edges%field(name) )
  end subroutine create_vector_field_in_edges

  subroutine mark_output(name, geom)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: geom
    type(Field_type) :: field
    integer :: jname
    field = geom%fields%field(name)
    call geom%output_fields%add_field(field)
  end subroutine mark_output

  function scalar_field(name, geom) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: geom
    type(Field_type) :: field
    real(kind=jprw), pointer :: array(:)
    field = geom%fields%field(name)
    array => field%data1()
  end function scalar_field

  function vector_field(name, geom) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: geom
    type(Field_type) :: field
    real(kind=jprw), pointer :: array(:,:)
    field = geom%fields%field(name)
    array => field%data2()
  end function vector_field

  subroutine synchronise_array_rank1(array, geom)
    implicit none
    real(kind=jprw), dimension(:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: geom
    call geom%nodes_comm%synchronise( array )
  end subroutine synchronise_array_rank1

  subroutine synchronise_array_rank2(array, geom)
    implicit none
    real(kind=jprw), dimension(:,:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: geom
    call geom%nodes_comm%synchronise( array )
  end subroutine synchronise_array_rank2

  subroutine gather_array_rank1(array_loc, array_full, geom)
    implicit none
    real(kind=jprw), dimension(:), intent(in) :: array_loc
    real(kind=jprw), dimension(:), allocatable, intent(out) :: array_full    
    type(DataStructure_type), intent(inout) :: geom
    call geom%nodes_comm%gather( array_loc, array_full )
  end subroutine gather_array_rank1

  subroutine gather_array_rank2(array_loc, array_full, geom)
    implicit none
    real(kind=jprw), dimension(:,:), intent(in) :: array_loc
    real(kind=jprw), dimension(:,:), allocatable, intent(out) :: array_full
    type(DataStructure_type), intent(inout) :: geom
    call geom%nodes_comm%gather( array_loc, array_full )
  end subroutine gather_array_rank2

end module datastruct_module
