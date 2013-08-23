! Grid module
! -----------------
! This module contains the Grid class which contains
! the nodes, elements, and element type
module grid_module
  ! Base classes, specifying the interface
  use elements_module,   only : Element, ShapeFunction
  
  ! All possible implementations for Element and Shapefunction respectively
  use lagrangep0_module, only : &
    & LagrangeP0_Quad2D,  LagrangeP0_Quad, &
    & LagrangeP0_Triag2D, LagrangeP0_Triag
  use lagrangep1_module, only : &
    & LagrangeP1_Quad2D,  LagrangeP1_Quad, &
    & LagrangeP1_Triag2D, LagrangeP1_Triag

  implicit none
  public
  
  type, public :: Field
    character(len=30) :: name
    real, dimension(:,:), allocatable :: array
    integer :: size
    integer :: cols
  contains
    procedure, pass :: init => Field__init
    procedure, pass :: destruct => Field__destruct
  end type Field
  
  type, public :: FieldPtr
    type(Field), pointer :: ptr
  end type FieldPtr
  
  type, public :: FunctionSpace
    class(Grid), pointer :: g
    class(ShapeFunction), allocatable :: sf
    character(len=30) :: name
    integer :: nb_elems
    integer :: nb_nodes
    integer, dimension(:,:), allocatable :: elements
    type(FieldPtr), dimension(:), allocatable :: fields
    
  contains
    procedure, pass :: init => FunctionSpace__init
    procedure, pass :: destruct => FunctionSpace__destruct
    procedure, pass :: add_field => FunctionSpace__add_field
    procedure, pass :: add_scalar_field => FunctionSpace__add_scalar_field
    procedure, pass :: add_vector_field => FunctionSpace__add_vector_field
  end type FunctionSpace
  
  type, public :: FunctionSpacePtr
    class(FunctionSpace), pointer :: ptr
  end type FunctionSpacePtr
  
  
  type, public, extends(FunctionSpace) :: ContinuousFunctionSpace
    
  contains
    procedure, pass :: init => ContinuousFunctionSpace__init
!     procedure, pass :: destruct => FunctionSpace__destruct
!     procedure, pass :: add_scalar_field => FunctionSpace__add_scalar_field
  end type ContinuousFunctionSpace
  
  type, public, extends(FunctionSpace) :: DiscontinuousFunctionSpace
  contains
    procedure, pass :: init => DiscontinuousFunctionSpace__init
!     procedure, pass :: destruct => FunctionSpace__destruct
!     procedure, pass :: add_scalar_field => FunctionSpace__add_scalar_field
  end type DiscontinuousFunctionSpace
  
  
  type, public :: Grid
    class(Element), allocatable :: elem
    integer :: dimension
    integer :: nb_elems
    integer :: nb_nodes
    real, dimension(:,:), allocatable    :: nodes
    integer, dimension(:,:), allocatable :: elements
    type(FunctionSpacePtr), dimension(:), allocatable :: function_spaces

  contains
    procedure, pass :: init => Grid__init
    procedure, pass :: destruct => Grid__destruct
    procedure, pass :: print => Grid__print
    procedure, pass :: add_discontinuous_function_space => Grid__add_discontinuous_function_space
    procedure, pass :: add_continuous_function_space => Grid__add_continuous_function_space
    procedure, pass :: elem_coords => Grid__elem_coords
  end type Grid
  
contains
  
  subroutine Field__init(self, name, function_space, cols)
    class(Field), intent(inout) :: self
    character(len=*), intent(in) :: name
    class(FunctionSpace), intent(in) :: function_space
    integer, intent(in) :: cols
    write(0,*) "Field::init(",name,",",function_space%name,")"  
    self%name = name
    self%cols = cols
    allocate(self%array(function_space%nb_nodes,self%cols))
    self%size = function_space%nb_nodes
    write(0,*) "Field::inited(",name,",",function_space%name,")"  
  end subroutine Field__init
  
  subroutine Field__destruct(self)
    class(Field), intent(inout) :: self
    write(0,*) "Field::destruct"
    if( allocated(self%array) ) deallocate(self%array)
  end subroutine Field__destruct
  
  subroutine FunctionSpace__init(self, name, shapefunction_type, grid_)
    class(FunctionSpace), intent(inout) :: self
    class(Grid), intent(in), target :: grid_
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "FunctionSpace::init(",shapefunction_type,")"
    self%name = name
    self%g => grid_
    ! Depending on the string element_type, the element will be allocated
    ! to be of a different type
    select case (shapefunction_type)
      case ("LagrangeP0")
        select case (grid_%elem%sf%shape)
          case (1) ! Quad
            allocate(LagrangeP0_Quad :: self%sf)    
          case (2) ! Triag
            allocate(LagrangeP0_Triag :: self%sf)
        end select
      case ("LagrangeP1")
        select case (grid_%elem%sf%shape)
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
    self%nb_elems = grid_%nb_elems
    allocate(self%fields(0))
  end subroutine FunctionSpace__init
  
  subroutine FunctionSpace__destruct(self)
    class(FunctionSpace), intent(inout) :: self
    integer :: f
    write(0,*) "FunctionSpace::destruct ",self%name
    if( allocated(self%sf) ) then
      call self%sf%destruct()
      deallocate(self%sf)
    endif
    if( allocated(self%elements) ) deallocate(self%elements)
    do f = 1,size(self%fields)
      if( associated(self%fields(f)%ptr) ) then
        call self%fields(f)%ptr%destruct()
        deallocate(self%fields(f)%ptr)
      endif
    end do
    deallocate(self%fields)
  end subroutine FunctionSpace__destruct
  
    function FunctionSpace__add_field(self,name,cols) result(new_field)
    class(FunctionSpace), intent(inout)   :: self
    character(len=*), intent(in) :: name
    integer, intent(in) :: cols
    type(Field), pointer :: new_field
    type(FieldPtr), allocatable :: tmp(:)
    allocate(new_field)
    call new_field%init(name,self,cols)
    call move_alloc(self%fields,tmp)
    allocate(self%fields(size(tmp)+1))
    self%fields(:size(tmp)) = tmp
    self%fields(size(self%fields))%ptr => new_field
    deallocate(tmp)
  end function FunctionSpace__add_field

  function FunctionSpace__add_scalar_field(self,name) result(new_field)
    class(FunctionSpace), intent(inout)   :: self
    character(len=*), intent(in) :: name
    type(Field), pointer :: new_field
    type(FieldPtr), allocatable :: tmp(:)
    allocate(new_field)
    call new_field%init(name,self,1)
    call move_alloc(self%fields,tmp)
    allocate(self%fields(size(tmp)+1))
    self%fields(:size(tmp)) = tmp
    self%fields(size(self%fields))%ptr => new_field
    deallocate(tmp)
  end function FunctionSpace__add_scalar_field
  
  function FunctionSpace__add_vector_field(self,name) result(new_field)
    class(FunctionSpace), intent(inout)   :: self
    character(len=*), intent(in) :: name
    type(Field), pointer :: new_field
    type(FieldPtr), allocatable :: tmp(:)
    allocate(new_field)
    call new_field%init(name,self,self%g%dimension)
    call move_alloc(self%fields,tmp)
    allocate(self%fields(size(tmp)+1))
    self%fields(:size(tmp)) = tmp
    self%fields(size(self%fields))%ptr => new_field
    deallocate(tmp)
  end function FunctionSpace__add_vector_field

  subroutine ContinuousFunctionSpace__init(self, name, shapefunction_type, grid_)
    class(ContinuousFunctionSpace), intent(inout) :: self
    class(Grid), intent(in), target :: grid_
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "ContinuousFunctionSpace::init(",shapefunction_type,")"
    call self%FunctionSpace%init(name,shapefunction_type,grid_)    
    self%nb_nodes = grid_%nb_nodes
    allocate(self%elements(self%nb_elems,self%sf%nb_nodes))
  end subroutine ContinuousFunctionSpace__init
  
  subroutine DiscontinuousFunctionSpace__init(self, name, shapefunction_type, grid_)
    class(DiscontinuousFunctionSpace), intent(inout) :: self
    class(Grid), intent(in), target :: grid_
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    write(0,*) "DiscontinuousFunctionSpace::init(",shapefunction_type,")"
    call self%FunctionSpace%init(name,shapefunction_type,grid_)    
    self%nb_nodes = self%nb_elems * self%sf%nb_nodes
    allocate(self%elements(self%nb_elems,self%sf%nb_nodes))
  end subroutine DiscontinuousFunctionSpace__init
  
  subroutine Grid__init(self, element_type)
    class(Grid), intent(inout) :: self
    character(len=*) :: element_type
    write(0,*) "Grid::init(",element_type,")"
    
    ! Depending on the string element_type, the element will be allocated
    ! to be of a different type
    select case (element_type)
      case ("LagrangeP1_Quad2D")
        allocate(LagrangeP1_Quad2D :: self%elem)
      case ("LagrangeP1_Triag2D")
        allocate(LagrangeP1_Triag2D :: self%elem)
      case default
        write(0,*) "Element type ",element_type," not supported"
        stop 1        
    end select
    call self%elem%init()
    self%dimension = self%elem%dimension
    allocate( FunctionSpacePtr :: self%function_spaces(0))
  end subroutine Grid__init
  
  subroutine Grid__destruct(self)
    class(Grid), intent(inout) :: self
    integer :: f
    write(0,*) "Grid::destruct"
    if( allocated(self%elem) ) then
      call self%elem%destruct()
      deallocate(self%elem)
    endif
    if( allocated(self%nodes) ) deallocate(self%nodes)
    if( allocated(self%elements) ) deallocate(self%elements)
    do f = 1,size(self%function_spaces)
      if( associated(self%function_spaces(f)%ptr) ) then
        call self%function_spaces(f)%ptr%destruct()
        deallocate(self%function_spaces(f)%ptr)
      endif
    end do
    deallocate(self%function_spaces)
  end subroutine Grid__destruct
  
  subroutine Grid__print(self)
    class(Grid), intent(inout) :: self
    integer :: n
    integer :: e
    write(0,*) "nodes = "
    do n = 1,self%nb_nodes
      write(0,*) "  ",self%nodes(n,:)
    end do
    write(0,*) "elems = "
    do e = 1,self%nb_elems
      write(0,*) "  ",self%elements(e,:)
    end do
  end subroutine Grid__print
  
  function Grid__add_continuous_function_space(self,name,shapefunction_type) result(new_function_space)
    class(Grid), intent(inout)   :: self
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    class(FunctionSpace), pointer :: new_function_space
    type(FunctionSpacePtr), allocatable :: tmp(:)
    integer :: f
    allocate( ContinuousFunctionSpace :: new_function_space )
    call new_function_space%init(name,shapefunction_type,self)
    call move_alloc(from=self%function_spaces,to=tmp)
    allocate(self%function_spaces(size(tmp)+1))
    do f=1,size(tmp)
      self%function_spaces(f)%ptr => tmp(f)%ptr
    end do
    self%function_spaces(size(tmp)+1)%ptr => new_function_space
  end function Grid__add_continuous_function_space

  function Grid__add_discontinuous_function_space(self,name,shapefunction_type) result(new_function_space)
    class(Grid), intent(inout)   :: self
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: shapefunction_type
    class(FunctionSpace), pointer :: new_function_space
    type(FunctionSpacePtr), allocatable :: tmp(:)
    integer :: f
    allocate( DiscontinuousFunctionSpace :: new_function_space )
    call new_function_space%init(name,shapefunction_type,self)
    call move_alloc(from=self%function_spaces,to=tmp)
    allocate(self%function_spaces(size(tmp)+1))
    do f=1,size(tmp)
      self%function_spaces(f)%ptr => tmp(f)%ptr
    end do
    self%function_spaces(size(tmp)+1)%ptr => new_function_space
  end function Grid__add_discontinuous_function_space

  subroutine Grid__elem_coords(self, elem_idx, elem_coords)
    class(Grid), intent(in) :: self
    integer, intent(in) :: elem_idx
    real, dimension(:,:), intent(inout) :: elem_coords

    integer n
    do n=1,self%elem%nb_nodes
      elem_coords(n,:) = self%nodes( self%elements(elem_idx,n), : )
    end do
  end subroutine Grid__elem_coords

end module grid_module