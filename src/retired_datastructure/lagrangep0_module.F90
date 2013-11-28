! LagrangeP0 module
! -----------------
! This module extends the elements_module module with concrete implementations
! for Lagrange P0 shape functions and element types
module lagrangep0_module
  use common_module, only: jprw, log_debug
  use elements_module, only: Element, Shapefunction_class
  implicit none
  public
  
  type, public, extends(ShapeFunction_class) :: LagrangeP0_Line
  contains
    procedure, public :: init        => LagrangeP0_Line__init
    procedure, public :: values      => LagrangeP0_Line__values
    procedure, public :: grad_values => LagrangeP0_Line__grad_values
  end type LagrangeP0_Line
  
  type, public, extends(Element) :: LagrangeP0_Line2D
  contains
    procedure, public :: init                 => LagrangeP0_Line2D__init
    procedure, public :: jacobian             => LagrangeP0_Line2D__jacobian
    procedure, public :: jacobian_determinant => LagrangeP0_Line2D__jacobian_determinant
  end type LagrangeP0_Line2D

  type, public, extends(ShapeFunction_class) :: LagrangeP0_Quad
  contains
    procedure, public :: init        => LagrangeP0_Quad__init
    procedure, public :: values      => LagrangeP0_Quad__values
    procedure, public :: grad_values => LagrangeP0_Quad__grad_values

  end type LagrangeP0_Quad
  
  type, public, extends(Element) :: LagrangeP0_Quad2D
  contains
    procedure, public :: init => LagrangeP0_Quad2D__init
  end type LagrangeP0_Quad2D
  
  
  type, public, extends(ShapeFunction_class) :: LagrangeP0_Triag
  contains
    procedure, public :: init        => LagrangeP0_Triag__init
    procedure, public :: values      => LagrangeP0_Triag__values
    procedure, public :: grad_values => LagrangeP0_Triag__grad_values


  end type LagrangeP0_Triag
  
  type, public, extends(Element) :: LagrangeP0_Triag2D
  contains
    procedure, public :: init => LagrangeP0_Triag2D__init
  end type LagrangeP0_Triag2D
  
contains

  subroutine LagrangeP0_Line__init(self)
    class(LagrangeP0_Line), intent(inout) :: self
    ! Call parent class init
    call self%ShapeFunction_class%init()
    self%shape = 0 ! Line
    self%dimensionality = 1
    self%nb_nodes = 1
    self%nb_sides = 2
    self%order = 0
    ! Deallocation will happen in base class ShapeFunction__destruct
    !allocate( self%local_coords(self%nb_nodes,self%dimensionality) )
    !self%local_coords(1,1) = 0.
    
    ! print something
    call log_debug("LagrangeP0_Line::init")
  end subroutine LagrangeP0_Line__init  

  subroutine LagrangeP0_Line__values( self, local_coord, values )
    class(LagrangeP0_Line), intent(in) :: self
    real(kind=jprw), dimension(:), intent(in) :: local_coord
    real(kind=jprw), dimension(:), intent(inout) :: values
    values(1) = 1.
    call log_debug("LagrangeP0_Line::values")
  end subroutine LagrangeP0_Line__values

  subroutine LagrangeP0_Line__grad_values( self, local_coord, grad_values )
    class(LagrangeP0_Line), intent(in) :: self
    real(kind=jprw), dimension(:), intent(in) :: local_coord
    real(kind=jprw), dimension(:,:), intent(inout) :: grad_values

    grad_values(1,1) = 0.
    grad_values(1,2) = 0.

    call log_debug( "LagrangeP0_Line::grad_values" )
  end subroutine LagrangeP0_Line__grad_values

  subroutine LagrangeP0_Line2D__init(self)
    class(LagrangeP0_Line2D), intent(inout) :: self

    ! Specify self%sf
    allocate(LagrangeP0_Line :: self%sf)
    call self%sf%init()

    ! Call parent class init
    call self%Element%init()    

    self%dimension = 2
    
    ! print something
    call log_debug( "LagrangeP0_Line2D::init" )
  end subroutine LagrangeP0_Line2D__init
    
    subroutine LagrangeP0_Line2D__jacobian( self, local_coord, elem_coords, jacobian )
    class(LagrangeP0_Line2D), intent(in)     :: self
    real(kind=jprw), dimension(:),        intent(in)    :: local_coord
    real(kind=jprw), dimension(:,:),      intent(in)    :: elem_coords
    real(kind=jprw), dimension(:,:),      intent(inout) :: jacobian
    call log_debug( "LagrangeP0_Line2D::jacobian()" )
  end subroutine LagrangeP0_Line2D__jacobian

  subroutine LagrangeP0_Line2D__jacobian_determinant( self, local_coord, elem_coords, jacobian_determinant )
    class(LagrangeP0_Line2D), intent(in)   :: self
    real(kind=jprw), dimension(:),       intent(in)   :: local_coord
    real(kind=jprw), dimension(:,:),     intent(in)   :: elem_coords
    real(kind=jprw), intent(out) :: jacobian_determinant
    jacobian_determinant = 0.
    write(0,*) "LagrangeP0_Line2D::jacobian_determinant()"
  end subroutine LagrangeP0_Line2D__jacobian_determinant


  subroutine LagrangeP0_Quad__init(self)
    class(LagrangeP0_Quad), intent(inout) :: self
    ! Call parent class init
    call self%ShapeFunction_class%init()
    self%shape = 1 ! Quad
    self%dimensionality = 2
    self%nb_nodes = 1
    self%nb_sides = 4
    self%order = 0
    ! Deallocation will happen in base class ShapeFunction__destruct
    allocate( self%local_coords(self%nb_nodes,self%dimensionality) )
    self%local_coords(1,:) = [0,0]

    ! print something
    call log_debug( "LagrangeP0_Quad::init" )
  end subroutine LagrangeP0_Quad__init

  subroutine LagrangeP0_Quad__values( self, local_coord, values )
    class(LagrangeP0_Quad), intent(in) :: self
    real(kind=jprw), dimension(:), intent(in) :: local_coord
    real(kind=jprw), dimension(:), intent(inout) :: values
    values(1)=1.
    write(0,*) "LagrangeP0_Quad::values"
  end subroutine LagrangeP0_Quad__values

  subroutine LagrangeP0_Quad__grad_values( self, local_coord, grad_values )
    class(LagrangeP0_Quad), intent(in) :: self
    real(kind=jprw), dimension(:), intent(in) :: local_coord
    real(kind=jprw), dimension(:,:), intent(inout) :: grad_values

    grad_values(1,1) = 0.
    grad_values(2,1) = 0.

    call log_debug( "LagrangeP0_Quad::grad_values" )
  end subroutine LagrangeP0_Quad__grad_values

  subroutine LagrangeP0_Quad2D__init(self)
    class(LagrangeP0_Quad2D), intent(inout) :: self

    ! Specify self%sf
    allocate(LagrangeP0_Quad :: self%sf)
    call self%sf%init()

    ! Call parent class init
    call self%Element%init()    

    self%dimension = 2
    
    ! print something
    call log_debug( "LagrangeP0_Quad2D::init" )
  end subroutine LagrangeP0_Quad2D__init
    
    
  subroutine LagrangeP0_Triag__init(self)
    class(LagrangeP0_Triag), intent(inout) :: self
    ! Call parent class init
    call self%ShapeFunction_class%init()

    self%shape = 2 ! Triag
    self%dimensionality = 2
    self%nb_nodes = 1
    self%nb_sides = 3
    self%order = 0
    ! Deallocation will happen in base class ShapeFunction__destruct
    allocate( self%local_coords(self%nb_nodes,self%dimensionality) )
    self%local_coords(1,:) = [1./3., 1./3.]

    ! print something
    call log_debug( "LagrangeP0_Triag::init" )
  end subroutine LagrangeP0_Triag__init  

  subroutine LagrangeP0_Triag__values( self, local_coord, values )
    class(LagrangeP0_Triag), intent(in) :: self
    real(kind=jprw), dimension(:), intent(in) :: local_coord
    real(kind=jprw), dimension(:), intent(inout) :: values
    values(1)=1.
    call log_debug( "LagrangeP0_Triag::values" )
  end subroutine LagrangeP0_Triag__values

  subroutine LagrangeP0_Triag__grad_values( self, local_coord, grad_values )
    class(LagrangeP0_Triag), intent(in) :: self
    real(kind=jprw), dimension(:), intent(in) :: local_coord
    real(kind=jprw), dimension(:,:), intent(inout) :: grad_values

    grad_values(1,1) = 0.
    grad_values(2,1) = 0.

    call log_debug( "LagrangeP0_Triag::grad_values" )
  end subroutine LagrangeP0_Triag__grad_values

  subroutine LagrangeP0_Triag2D__init(self)
    class(LagrangeP0_Triag2D), intent(inout) :: self

    ! Specify self%sf
    allocate(LagrangeP0_Triag :: self%sf)
    call self%sf%init()

    ! Call parent class init
    call self%Element%init()    

    self%dimension = 2
    
    ! print something
    call log_debug( "LagrangeP0_Triag2D::init" )
  end subroutine LagrangeP0_Triag2D__init
  
end module lagrangep0_module
