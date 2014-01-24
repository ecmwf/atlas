module datastruct_module

  use common_module, only: jprw
  use datastruct, only : Mesh_type, FieldSet_type, Field_type, HaloExchange_type, &
      & new_FunctionSpace, new_PrismaticFunctionSpace, FunctionSpace_type, new_Mesh, new_FieldSet
  use parallel_module

  implicit none
  private
  public :: create_mesh, create_mesh_3d
  public :: create_field_in_nodes_2d
  public :: create_field_in_edges_2d
  public :: create_field_in_nodes_3d
  public :: create_field_in_edges_3d
  public :: scalar_field_2d, vector_field_2d
  public :: scalar_field_3d, vector_field_3d
  public :: mark_output
  public :: gather

  public :: halo_exchange, halo_exchange_2d, halo_exchange_3d

  interface halo_exchange
    module procedure halo_exchange_2d_real_rank1
    module procedure halo_exchange_2d_real_rank2
  end interface halo_exchange

  interface halo_exchange_2d
    module procedure halo_exchange_2d_real_rank1
    module procedure halo_exchange_2d_real_rank2
  end interface halo_exchange_2d

  interface halo_exchange_3d
    module procedure halo_exchange_3d_real_rank2
    module procedure halo_exchange_3d_real_rank3
  end interface halo_exchange_3d



  interface  gather
    module procedure gather_array_rank1
    module procedure gather_array_rank2
    module procedure gather_array_rank3
  end interface gather
!  type, public, extends(Mesh_type) :: Mesh_class
!  end type Mesh_class

  type, public :: DataStructure_type
  private
  ! PRIVATE part
  ! ------------
    type(Mesh_type)     :: mesh
    type(FieldSet_type) :: fields
    type(FunctionSpace_type) :: functionspace_nodes_2d
    type(FunctionSpace_type) :: functionspace_edges_2d
    type(FunctionSpace_type) :: functionspace_nodes_3d
    type(FunctionSpace_type) :: functionspace_edges_3d
    type(Comm_type), public :: nodes_comm

    type(FieldSet_type), public :: output_fields

    real(kind=jprw), public :: time

  ! PUBLIC part
  ! -----------
    integer,                       public :: glb_nb_nodes=0
    integer,                       public :: glb_nb_edges=0
    integer,                       public :: nb_nodes=0
    integer,                       public :: nb_edges=0
    integer,                       public :: nb_levels=0
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
  
  subroutine create_mesh(nb_nodes, nb_edges, proc, glb_idx, dstruct)
    implicit none
    integer, intent(in) :: nb_nodes
    integer, intent(in) :: nb_edges
    integer, intent(in) :: proc(:)
    integer, intent(in) :: glb_idx(:)
    type(DataStructure_type), intent(inout), target :: dstruct

    dstruct%mesh = new_Mesh()
    dstruct%fields = new_FieldSet("fields")
    dstruct%output_fields = new_FieldSet("output_fields")

    call dstruct%mesh%add_function_space( new_FunctionSpace("nodes_2d", "P1-C", nb_nodes) )
    dstruct%functionspace_nodes_2d = dstruct%mesh%function_space("nodes_2d")
    dstruct%nb_nodes = nb_nodes

    call dstruct%mesh%add_function_space( new_FunctionSpace("edges_2d", "P0-D", nb_edges) )
    dstruct%functionspace_edges_2d = dstruct%mesh%function_space("edges_2d")
    dstruct%nb_edges = nb_edges

    allocate( dstruct%edges(dstruct%nb_edges,2) )

    call create_field_in_nodes_2d("coordinates",2,dstruct)
    call create_field_in_nodes_2d("dual_volumes",1,dstruct)
    call create_field_in_edges_2d("dual_normals",2,dstruct)

    allocate( dstruct%nodes_proc( nb_nodes ) ) 
    allocate( dstruct%nodes_glb_idx( nb_nodes ) ) 
    dstruct%nodes_proc(:) = proc(:)
    dstruct%nodes_glb_idx(:) = glb_idx(:)
    call dstruct%nodes_comm%setup( proc, glb_idx )

    call dstruct%functionspace_nodes_2d%parallelise(proc,glb_idx)

    allocate( dstruct%edges_proc( nb_edges ) ) 
    allocate( dstruct%edges_glb_idx( nb_edges ) ) 

  end subroutine create_mesh

  subroutine create_mesh_3d(nb_levels, nb_nodes, nb_edges, proc, glb_idx, dstruct)
    implicit none
    integer, intent(in) :: nb_levels
    integer, intent(in) :: nb_nodes
    integer, intent(in) :: nb_edges
    integer, intent(in) :: proc(:)
    integer, intent(in) :: glb_idx(:)
    type(DataStructure_type), intent(inout), target :: dstruct

    dstruct%mesh = new_Mesh()
    dstruct%fields = new_FieldSet("fields")
    dstruct%output_fields = new_FieldSet("output_fields")

    dstruct%nb_levels = nb_levels

    dstruct%nb_nodes = nb_nodes
    call dstruct%mesh%add_function_space( new_FunctionSpace("nodes_2d", "P1-C", nb_nodes) )
    dstruct%functionspace_nodes_2d = dstruct%mesh%function_space("nodes_2d")

    call dstruct%mesh%add_function_space( new_PrismaticFunctionSpace("nodes_3d", "P1-C", nb_levels, nb_nodes) )
    dstruct%functionspace_nodes_3d = dstruct%mesh%function_space("nodes_3d")

    call dstruct%mesh%add_function_space( new_FunctionSpace("edges_2d", "P0-D", nb_edges) )
    dstruct%functionspace_edges_2d = dstruct%mesh%function_space("edges_2d")

    call dstruct%mesh%add_function_space( new_PrismaticFunctionSpace("edges_3d", "P0-D", nb_levels, nb_edges) )
    dstruct%functionspace_edges_3d = dstruct%mesh%function_space("edges_3d")
    dstruct%nb_edges = nb_edges

    allocate( dstruct%edges(dstruct%nb_edges,2) )

    call create_field_in_nodes_2d("coordinates",2,dstruct)
    call create_field_in_nodes_3d("dual_volumes",1,dstruct)
    call create_field_in_edges_2d("dual_normals",2,dstruct)

    allocate( dstruct%nodes_proc( nb_nodes ) ) 
    allocate( dstruct%nodes_glb_idx( nb_nodes ) ) 
    dstruct%nodes_proc(:) = proc(:)
    dstruct%nodes_glb_idx(:) = glb_idx(:)
    call dstruct%nodes_comm%setup( proc, glb_idx )

    allocate( dstruct%edges_proc( nb_edges ) ) 
    allocate( dstruct%edges_glb_idx( nb_edges ) ) 

    call dstruct%functionspace_nodes_2d%parallelise(proc,glb_idx)
    call dstruct%functionspace_nodes_3d%parallelise(proc,glb_idx)

  end subroutine create_mesh_3d

  subroutine create_field_in_nodes_3d(name, nb_vars, dstruct)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: nb_vars
    type(DataStructure_type), intent(inout) :: dstruct
    type(FunctionSpace_type) :: function_space
    function_space = dstruct%functionspace_nodes_3d
    call function_space%create_field(name,nb_vars)
    call dstruct%fields%add_field( function_space%field(name) )
  end subroutine create_field_in_nodes_3d

  subroutine create_field_in_edges_3d(name, nb_vars, dstruct)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: nb_vars
    type(DataStructure_type), intent(inout) :: dstruct
    type(FunctionSpace_type) :: function_space
    function_space = dstruct%functionspace_edges_3d
    call function_space%create_field(name,nb_vars)
    call dstruct%fields%add_field( function_space%field(name) )
  end subroutine create_field_in_edges_3d

  subroutine create_field_in_nodes_2d(name, nb_vars, dstruct)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: nb_vars
    type(DataStructure_type), intent(inout) :: dstruct
    type(FunctionSpace_type) :: function_space
    function_space = dstruct%functionspace_nodes_2d
    call function_space%create_field(name,nb_vars)
    call dstruct%fields%add_field( function_space%field(name) )
  end subroutine create_field_in_nodes_2d

  subroutine create_field_in_edges_2d(name, nb_vars, dstruct)
    implicit none
    character(len=*), intent(in) :: name
    integer, intent(in) :: nb_vars
    type(DataStructure_type), intent(inout) :: dstruct
    type(FunctionSpace_type) :: function_space
    function_space = dstruct%functionspace_edges_2d
    call function_space%create_field(name,nb_vars)
    call dstruct%fields%add_field( function_space%field(name) )
  end subroutine create_field_in_edges_2d

  subroutine mark_output(name, dstruct)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: dstruct
    type(Field_type) :: field
    field = dstruct%fields%field(name)
    call dstruct%output_fields%add_field(field)
  end subroutine mark_output

  function scalar_field_2d(name, dstruct) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: dstruct
    type(Field_type) :: field
    real(kind=jprw), pointer :: array(:)
    field = dstruct%fields%field(name)
    array => field%data1()
  end function scalar_field_2d

  function vector_field_2d(name, dstruct) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: dstruct
    type(Field_type) :: field
    real(kind=jprw), pointer :: array(:,:)
    field = dstruct%fields%field(name)
    array => field%data2()
  end function vector_field_2d

  function scalar_field_3d(name, dstruct) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: dstruct
    type(Field_type) :: field
    real(kind=jprw), pointer :: array(:,:)
    field = dstruct%fields%field(name)
    array => field%data2()
  end function scalar_field_3d

  function vector_field_3d(name, dstruct) result(array)
    implicit none
    character(len=*), intent(in) :: name
    type(DataStructure_type), intent(in) :: dstruct
    type(Field_type) :: field
    real(kind=jprw), pointer :: array(:,:,:)
    field = dstruct%fields%field(name)
    array => field%data3()
  end function vector_field_3d

  subroutine halo_exchange_2d_real_rank1(array, dstruct)
    implicit none
    real(kind=jprw), dimension(:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%functionspace_nodes_2d%halo_exchange( array )
  end subroutine halo_exchange_2d_real_rank1

   subroutine halo_exchange_2d_real_rank2(array, dstruct)
    implicit none
    real(kind=jprw), dimension(:,:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%functionspace_nodes_2d%halo_exchange( array )
  end subroutine halo_exchange_2d_real_rank2

   subroutine halo_exchange_3d_real_rank2(array, dstruct)
    implicit none
    real(kind=jprw), dimension(:,:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%functionspace_nodes_3d%halo_exchange( array )
  end subroutine halo_exchange_3d_real_rank2

  subroutine halo_exchange_3d_real_rank3(array, dstruct)
    implicit none
    real(kind=jprw), dimension(:,:,:), intent(inout) :: array
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%functionspace_nodes_3d%halo_exchange( array )
  end subroutine halo_exchange_3d_real_rank3

  subroutine gather_array_rank1(array_loc, array_full, dstruct)
    implicit none
    real(kind=jprw), dimension(:), intent(in) :: array_loc
    real(kind=jprw), dimension(:), allocatable, intent(out) :: array_full    
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%nodes_comm%gather( array_loc, array_full )
  end subroutine gather_array_rank1

  subroutine gather_array_rank2(array_loc, array_full, dstruct)
    implicit none
    real(kind=jprw), dimension(:,:), intent(in) :: array_loc
    real(kind=jprw), dimension(:,:), allocatable, intent(out) :: array_full
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%nodes_comm%gather( array_loc, array_full )
  end subroutine gather_array_rank2

  subroutine gather_array_rank3(array_loc, array_full, dstruct)
    implicit none
    real(kind=jprw), dimension(:,:,:), intent(in) :: array_loc
    real(kind=jprw), dimension(:,:,:), allocatable, intent(out) :: array_full
    type(DataStructure_type), intent(inout) :: dstruct
    call dstruct%nodes_comm%gather( array_loc, array_full )
  end subroutine gather_array_rank3

end module datastruct_module
