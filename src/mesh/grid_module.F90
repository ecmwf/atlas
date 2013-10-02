! Grid module
! -----------------
! This module contains the Grid class which contains
! the nodes, elements, and element type
module grid_module
  ! Base classes, specifying the interface
  use common_module
  use elements_module,   only : Element, ShapeFunction_class

  ! All possible implementations for Element and Shapefunction respectively
  use lagrangep0_module, only : &
    & LagrangeP0_Line2D,  LagrangeP0_Line, &
    & LagrangeP0_Quad2D,  LagrangeP0_Quad, &
    & LagrangeP0_Triag2D, LagrangeP0_Triag
  use lagrangep1_module, only : &
    & LagrangeP1_Line2D,  LagrangeP1_Line, &
    & LagrangeP1_Quad2D,  LagrangeP1_Quad, &
    & LagrangeP1_Triag2D, LagrangeP1_Triag

  implicit none
  public
  
  public new_FaceFunctionSpace
  public new_DiscontinuousFunctionSpace


  type, public :: Field_class
    character(len=30) :: name
    real(kind=jprb), dimension(:,:), allocatable :: array
    integer :: rows
    integer :: cols
    class(FunctionSpace_class), pointer :: function_space
  contains
    procedure, pass :: init => Field__init
    procedure, pass :: destruct => Field__destruct
  end type Field_class
  
  type, public :: FieldPtr
    type(Field_class), pointer :: ptr
  end type FieldPtr
  
  type, public :: FunctionSpace_class
    class(Grid_class), pointer :: grid
    class(ShapeFunction_class), allocatable :: sf
    character(len=30) :: name
    integer :: nb_elems = 0
    integer :: nb_nodes = 0
    integer :: nb_fields = 0
    integer, dimension(:,:), allocatable :: elements
    type(Field_class), dimension(:), pointer :: fields
    
  contains
    procedure, pass :: init => FunctionSpace__init
    procedure, pass :: destruct => FunctionSpace__destruct
    procedure, pass :: add_scalar_field => FunctionSpace__add_scalar_field
    procedure, pass :: add_vector_field => FunctionSpace__add_vector_field
    procedure, pass :: add_array_field => FunctionSpace__add_array_field
    procedure, pass :: field => FunctionSpace__field


  end type FunctionSpace_class
  
  type, public :: FunctionSpacePtr
    class(FunctionSpace_class), pointer :: ptr
  end type FunctionSpacePtr
  
  
  type, public, extends(FunctionSpace_class) :: ContinuousFunctionSpace
    
  contains
    procedure, pass :: init => ContinuousFunctionSpace__init
  end type ContinuousFunctionSpace
  
  type, public, extends(FunctionSpace_class) :: DiscontinuousFunctionSpace
  contains
    procedure, pass :: init => DiscontinuousFunctionSpace__init
  end type DiscontinuousFunctionSpace
  
  type, public, extends(FunctionSpace_class) :: FaceFunctionSpace
  contains
    procedure, pass :: init => FaceFunctionSpace__init
  end type FaceFunctionSpace


  type, public :: Grid_class
    class(Element), allocatable :: cell
    class(Element), allocatable :: face
    !type(ContinuousFunctionSpace) :: nodes
    integer :: dimension
    integer :: nb_elems
    integer :: nb_nodes
    integer :: nb_faces
    integer :: nb_pole_faces
    integer :: nb_internal_faces
    integer :: nb_periodic_nodes
    real(kind=jprb),  dimension(:,:), allocatable :: nodes_coordinates
    integer, dimension(:,:), allocatable :: cells
    integer, dimension(:,:), allocatable :: faces
    integer, dimension(:), allocatable :: internal_faces
    integer, dimension(:), allocatable :: pole_faces
    integer, dimension(:,:), allocatable :: periodic_nodes

    type(FunctionSpacePtr), dimension(:), allocatable :: function_spaces

  contains
    procedure, pass :: init => Grid__init
    procedure, pass :: destruct => Grid__destruct
    procedure, pass :: print => Grid__print
    procedure, pass :: cell_coords => Grid__cell_coords
    procedure, pass :: function_space => Grid__function_space
    procedure, pass :: add_function_space  => Grid__add_function_space

  end type Grid_class
  
  type, public :: State_class
    character(len=30) :: name
    real(kind=jprb)              :: time
    integer :: nb_fields = 0
    type(FieldPtr), dimension(:), allocatable :: fields
  contains
    procedure, pass :: init => State__init
    procedure, pass :: destruct => State__destruct
    procedure, pass :: add_field => State__add_field
    procedure, pass :: field => State__field
    procedure, pass :: has_field => State__has_field
  end type State_class
  
contains
  
! ------------------------------------------------------------------------------------
!                                  Field subroutines
!-------------------------------------------------------------------------------------

  subroutine Field__init(self, name, function_space, cols)
    class(Field_class), intent(inout) :: self
    character(len=*), intent(in) :: name
    class(FunctionSpace_class), intent(in), target :: function_space
    integer, intent(in) :: cols
    write(0,*) "Field::init(",name,",",function_space%name,"),"
    self%name = name
    self%function_space => function_space
    self%cols = cols
    allocate(self%array(self%function_space%nb_nodes,self%cols))
    self%rows = self%function_space%nb_nodes
  end subroutine Field__init
  
  subroutine Field__destruct(self)
    class(Field_class), intent(inout) :: self
    write(0,*) "Field::destruct"
    if( allocated(self%array) ) deallocate(self%array)
  end subroutine Field__destruct
  

! ------------------------------------------------------------------------------------
!                                 FunctionSpace subroutines
!-------------------------------------------------------------------------------------

  subroutine FunctionSpace__init(self, name, shapefunction_type, grid)
    class(FunctionSpace_class), intent(inout) :: self
    class(Grid_class), intent(in), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "FunctionSpace::init(",shapefunction_type,")"
    self%name = name
    self%grid => grid
    ! Depending on the string element_type, the element will be allocated
    ! to be of a different type
    select case (shapefunction_type)
      case ("LagrangeP0")
        select case (self%grid%cell%sf%shape)
          case (1) ! Quad
            allocate(LagrangeP0_Quad :: self%sf)    
          case (2) ! Triag
            allocate(LagrangeP0_Triag :: self%sf)
        end select
      case ("LagrangeP1")
        select case (self%grid%cell%sf%shape)
          case (1) ! Quad
            allocate(LagrangeP1_Quad :: self%sf)    
          case (2) ! Triag
            allocate(LagrangeP1_Triag :: self%sf)
        end select
      case default
        write(0,*) "ShapeFunction type ",shapefunction_type," not supported"
        stop 1        
    end select
    call self%sf%init()
    self%nb_elems = self%grid%nb_elems
    allocate(self%fields(100))
  end subroutine FunctionSpace__init
  
  subroutine FunctionSpace__destruct(self)
    class(FunctionSpace_class), intent(inout) :: self
    integer :: f
    write(0,*) "FunctionSpace::destruct ",self%name
    if( allocated(self%sf) ) then
      call self%sf%destruct()
      deallocate(self%sf)
    endif
    if( allocated(self%elements) ) deallocate(self%elements)
  end subroutine FunctionSpace__destruct
  
  function FunctionSpace__add_array_field(self,name,cols) result(new_field)
    class(FunctionSpace_class), intent(inout)   :: self
    character(len=*), intent(in) :: name
    integer, intent(in) :: cols
    class(Field_class), pointer :: new_field
    self%nb_fields = self%nb_fields+1
    new_field => self%fields(self%nb_fields)
    call new_field%init(name,self,cols)
  end function FunctionSpace__add_array_field


  function FunctionSpace__add_scalar_field(self,name) result(new_field)
    class(FunctionSpace_class), intent(inout)   :: self
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: new_field
    new_field => self%add_array_field(name,1)
  end function FunctionSpace__add_scalar_field
  
  function FunctionSpace__add_vector_field(self,name) result(new_field)
    class(FunctionSpace_class), intent(inout)   :: self
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: new_field
    new_field => self%add_array_field(name,self%grid%dimension)
  end function FunctionSpace__add_vector_field

  function FunctionSpace__field(self, name) result(field_)
    class(FunctionSpace_class), intent(in) :: self
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: field_
    integer :: f
    do f=1,self%nb_fields
      field_ => self%fields(f)
      if( field_%name == name) then
        return
      end if
    end do
    write(0,*) 'call abort '
    call abort
  end function FunctionSpace__field

  subroutine ContinuousFunctionSpace__init(self, name, shapefunction_type, grid)
    class(ContinuousFunctionSpace), intent(inout) :: self
    class(Grid_class), intent(in), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "ContinuousFunctionSpace::init(",shapefunction_type,")"
    call self%FunctionSpace_class%init(name,shapefunction_type,grid)
    self%nb_nodes = grid%nb_nodes
    allocate(self%elements(self%nb_elems,self%sf%nb_nodes))
  end subroutine ContinuousFunctionSpace__init
  
  subroutine DiscontinuousFunctionSpace__init(self, name, shapefunction_type, grid)
    class(DiscontinuousFunctionSpace), intent(inout) :: self
    class(Grid_class), intent(in), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "DiscontinuousFunctionSpace::init(",shapefunction_type,")"
    call self%FunctionSpace_class%init(name,shapefunction_type,grid)    
    self%nb_nodes = self%nb_elems * self%sf%nb_nodes
    allocate(self%elements(self%nb_elems,self%sf%nb_nodes))
  end subroutine DiscontinuousFunctionSpace__init


  subroutine FaceFunctionSpace__init(self, name, shapefunction_type, grid)
    class(FaceFunctionSpace), intent(inout) :: self
    class(Grid_class), intent(in), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "FaceFunctionSpace::init(",shapefunction_type,")"
    self%name = name
    self%grid => grid
    self%nb_elems = grid%nb_faces
    ! Depending on the string element_type, the element will be allocated
    ! to be of a different type
    select case (shapefunction_type)
      case ("LagrangeP0")
        self%nb_nodes = self%nb_elems
        select case (grid%face%sf%shape)
          case (0) ! Line
            allocate(LagrangeP0_Line :: self%sf)
          case (1) ! Quad
            allocate(LagrangeP0_Quad :: self%sf)    
          case (2) ! Triag
            allocate(LagrangeP0_Triag :: self%sf)
        end select
      case ("LagrangeP1")
        self%nb_nodes = grid%nb_nodes
        select case (grid%face%sf%shape)
          case (0) ! Line
            allocate(LagrangeP1_Line :: self%sf)
          case (1) ! Quad
            allocate(LagrangeP1_Quad :: self%sf)    
          case (2) ! Triag
            allocate(LagrangeP1_Triag :: self%sf)
        end select
      case default
        write(0,*) "ShapeFunction type ",shapefunction_type," not supported"
        stop 1        
    end select
    call self%sf%init()
    allocate(self%fields(100))
    allocate(self%elements(self%nb_elems,self%sf%nb_nodes))
  end subroutine FaceFunctionSpace__init

  function new_FaceFunctionSpace(name, shapefunction_type, grid) result(function_space)
    class(Grid_class), intent(inout), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    class(FunctionSpace_class), pointer :: function_space
    allocate( FaceFunctionSpace :: function_space )
    call function_space%init(name, shapefunction_type, grid)
    write(0,*) function_space%name
  end function new_FaceFunctionSpace

  function new_ContinuousFunctionSpace(name, shapefunction_type, grid) result(function_space)
    class(Grid_class), intent(inout), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    class(FunctionSpace_class), pointer :: function_space
    allocate( ContinuousFunctionSpace :: function_space )
    call function_space%init(name, shapefunction_type, grid)
    write(0,*) function_space%name
  end function new_ContinuousFunctionSpace

  function new_DiscontinuousFunctionSpace(name, shapefunction_type, grid) result(function_space)
    class(Grid_class), intent(inout), target :: grid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    class(FunctionSpace_class), pointer :: function_space
    allocate( DiscontinuousFunctionSpace :: function_space )
    call function_space%init(name, shapefunction_type, grid)
    write(0,*) function_space%name
  end function new_DiscontinuousFunctionSpace

  function new_ScalarField(name, function_space) result(field)
    class(FunctionSpace_class), intent(in), target :: function_space
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: field
    allocate(field)
    call field%init(name,function_space,1)
  end function new_ScalarField

  function new_VectorField(name, function_space) result(field)
    class(FunctionSpace_class), intent(in), target :: function_space
    character(len=*), intent(in) :: name
    class(Field_class), pointer :: field
    allocate(field)
    call field%init(name,function_space, function_space%grid%dimension)
  end function new_VectorField

! ------------------------------------------------------------------------------------
!                                   Grid subroutines
!-------------------------------------------------------------------------------------
  subroutine Grid__init(self, element_type)
    class(Grid_class), intent(inout) :: self
    character(len=*) :: element_type
    write(0,*) "Grid::init(",element_type,")"
    
    ! Depending on the string element_type, the element will be allocated
    ! to be of a different type
    select case (element_type)
      case ("LagrangeP1_Quad2D")
        allocate(LagrangeP1_Quad2D :: self%cell)
        allocate(LagrangeP1_Line2D :: self%face)
      case ("LagrangeP1_Triag2D")
        allocate(LagrangeP1_Triag2D :: self%cell)
        allocate(LagrangeP1_Line2D  :: self%face)
      case default
        write(0,*) "Element type ",element_type," not supported"
        stop 1        
    end select
    call self%cell%init()
    call self%face%init()
    self%dimension = self%cell%dimension
    allocate(self%function_spaces(0))
  end subroutine Grid__init
  
  subroutine Grid__destruct(self)
    class(Grid_class), intent(inout) :: self
    integer :: f
    write(0,*) "Grid::destruct"
    if( allocated(self%cell) ) then
      call self%cell%destruct()
      deallocate(self%cell)
    endif
    if( allocated(self%face) ) then
      call self%face%destruct()
      deallocate(self%face)
    endif
    if( allocated(self%nodes_coordinates) ) deallocate(self%nodes_coordinates)
    if( allocated(self%cells) ) deallocate(self%cells)
    if( allocated(self%faces) ) deallocate(self%faces)

  end subroutine Grid__destruct
  
  subroutine Grid__print(self)
    class(Grid_class), intent(inout) :: self
    integer :: n
    integer :: e
    write(0,*) "nodes = "
    do n = 1,self%nb_nodes
      write(0,*) "  ",self%nodes_coordinates(n,:)
    end do
    write(0,*) "elems = "
    do e = 1,self%nb_elems
      write(0,*) "  ",self%cells(e,:)
    end do
  end subroutine Grid__print
  
  subroutine Grid__cell_coords(self, elem_idx, cell_coords)
    class(Grid_class), intent(in) :: self
    integer, intent(in) :: elem_idx
    real(kind=jprb), dimension(:,:), intent(inout) :: cell_coords

    integer n
    do n=1,self%cell%nb_nodes
      cell_coords(n,:) = self%nodes_coordinates( self%cells(elem_idx,n), : )
    end do
  end subroutine Grid__cell_coords

  subroutine Grid__add_function_space(self,function_space)
    class(Grid_class), intent(inout) :: self
    class(FunctionSpace_class), target, intent(inout) :: function_space
    type(FunctionSpacePtr), allocatable :: tmp(:)
    integer :: f
    write(0,*) "adding ",function_space%name
    call move_alloc(from=self%function_spaces,to=tmp)
    allocate(self%function_spaces(size(tmp)+1))
    do f=1,size(tmp)
      self%function_spaces(f)%ptr => tmp(f)%ptr
    end do
    self%function_spaces(size(tmp)+1)%ptr => function_space
  end subroutine Grid__add_function_space


  function Grid__function_space(self, name) result(function_space)
    class(Grid_class), intent(in) :: self
    character(len=*), intent(in) :: name
    class(FunctionSpace_class), pointer :: function_space
    integer :: f
    write(0,*) "shape function_spaces = " , size(self%function_spaces)
    do f=1,size(self%function_spaces)
      function_space => self%function_spaces(f)%ptr
      if( function_space%name == name) then
        return
      end if
    end do
    write(0,*) "FunctionSpace ", trim(name), " not found"
    call abort
  end function Grid__function_space

! ------------------------------------------------------------------------------------
!                                  State subroutines
!-------------------------------------------------------------------------------------

  subroutine State__init(self, name)
    class(State_class), intent(inout) :: self
    character(len=*), intent(in) :: name
    write(0,*) "State::init(",name,")"  
    self%name = name  
    allocate(self%fields(100))
  end subroutine State__init
  
  subroutine State__destruct(self)
    class(State_class), intent(inout) :: self
    write(0,*) "State::destruct"
    deallocate(self%fields)
  end subroutine State__destruct
  
  subroutine State__add_field(self,name,function_space,cols)
    class(State_class), intent(inout) :: self
    character(len=*), intent(in) :: name
    class(FunctionSpace_class), intent(in), target ::function_space
    integer, intent(in) :: cols
    self%nb_fields = self%nb_fields + 1
    allocate( self%fields( self%nb_fields )%ptr )
    call self%fields( self%nb_fields )%ptr%init(name,function_space,cols)
  end subroutine State__add_field

  function State__field(self, name) result(field)
    class(State_class), intent(in) :: self
    character(len=*), intent(in) :: name
    type(Field_class), pointer :: field
    integer :: f
    do f=1,self%nb_fields
      field => self%fields(f)%ptr
      if( field%name == name) then
        return
      end if
    end do
    write(0,*) 'Could not find field ',name, ' in ',self%name
    call abort
  end function State__field

  logical function State__has_field(self, name)
    class(State_class), intent(in) :: self
    character(len=*), intent(in) :: name
    type(Field_class), pointer :: field
    logical :: has_field
    integer :: f
    do f=1,self%nb_fields
      field => self%fields(f)%ptr
      if( field%name == name) then
        State__has_field = .True.
        return
      end if
    end do
    State__has_field = .False.
  end function State__has_field

end module grid_module
