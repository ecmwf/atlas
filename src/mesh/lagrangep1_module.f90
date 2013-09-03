! LagrangeP1 module
! -----------------
! This module extends the elements_module module with concrete implementations
! for Lagrange P1 shape functions and element types
module lagrangep1_module
  use elements_module
  implicit none
  public
  
  type, public, extends(ShapeFunction_class) :: LagrangeP1_Line
  contains
    procedure, public :: init        => LagrangeP1_Line__init
    procedure, public :: values      => LagrangeP1_Line__values
    procedure, public :: grad_values => LagrangeP1_Line__grad_values
  end type LagrangeP1_Line
  
  type, public, extends(Element) :: LagrangeP1_Line2D
  contains
    procedure, public :: init                 => LagrangeP1_Line2D__init
    procedure, public :: jacobian             => LagrangeP1_Line2D__jacobian
    procedure, public :: jacobian_determinant => LagrangeP1_Line2D__jacobian_determinant
  end type LagrangeP1_Line2D

  type, public, extends(ShapeFunction_class) :: LagrangeP1_Quad
  contains
    procedure, public :: init        => LagrangeP1_Quad__init
    procedure, public :: values      => LagrangeP1_Quad__values
    procedure, public :: grad_values => LagrangeP1_Quad__grad_values
  end type LagrangeP1_Quad
  
  type, public, extends(Element) :: LagrangeP1_Quad2D
  contains
    procedure, public :: init             => LagrangeP1_Quad2D__init
    procedure, public :: jacobian         => LagrangeP1_Quad2D__jacobian
    procedure, public :: jacobian_inverse => LagrangeP1_Quad2D__jacobian_inverse
  end type LagrangeP1_Quad2D
  
  
  type, public, extends(ShapeFunction_class) :: LagrangeP1_Triag
  contains
    procedure, public :: init        => LagrangeP1_Triag__init
    procedure, public :: values      => LagrangeP1_Triag__values
    procedure, public :: grad_values => LagrangeP1_Triag__grad_values 
  end type LagrangeP1_Triag
  
  type, public, extends(Element) :: LagrangeP1_Triag2D
  contains
    procedure, public :: init             => LagrangeP1_Triag2D__init
    procedure, public :: jacobian         => LagrangeP1_Triag2D__jacobian
    procedure, public :: jacobian_inverse => LagrangeP1_Triag2D__jacobian_inverse

  end type LagrangeP1_Triag2D
  
contains

  subroutine LagrangeP1_Line__init(self)
    class(LagrangeP1_Line), intent(inout) :: self
    ! Call parent class init
    call self%ShapeFunction_class%init()
    self%shape = 0 ! Line
    self%dimensionality = 1
    self%nb_nodes = 2
    self%nb_sides = 2
    self%order = 1
    ! Deallocation will happen in base class ShapeFunction__destruct
    allocate( self%local_coords(self%nb_nodes,self%dimensionality) )
    self%local_coords(1,1) = -1.
    self%local_coords(2,1) = +1.

    ! print something
    write(0,*) "LagrangeP1_Line::init"
  end subroutine LagrangeP1_Line__init  

  subroutine LagrangeP1_Line__values( self, local_coord, values )
    class(LagrangeP1_Line), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:), intent(inout) :: values
    values(1) = 0.5 * (1.-local_coord(1))
    values(2) = 0.5 * (1.+local_coord(1))
    write(0,*) "LagrangeP1_Line::values"
  end subroutine LagrangeP1_Line__values

  subroutine LagrangeP1_Line__grad_values( self, local_coord, grad_values )
    class(LagrangeP1_Line), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:,:), intent(inout) :: grad_values

    grad_values(1,1) = -0.5
    grad_values(1,2) = +0.5

    write(0,*) "LagrangeP1_Line::grad_values"
  end subroutine LagrangeP1_Line__grad_values

    subroutine LagrangeP1_Line2D__init(self)
    class(LagrangeP1_Line2D), intent(inout) :: self

    ! Specify self%sf
    allocate(LagrangeP1_Line :: self%sf)
    call self%sf%init()

    ! Call parent class init
    call self%Element%init()    

    self%dimension = 2
    
    ! print something
    write(0,*) "LagrangeP1_Line2D::init"
  end subroutine LagrangeP1_Line2D__init
    
    subroutine LagrangeP1_Line2D__jacobian( self, local_coord, elem_coords, jacobian )
    class(LagrangeP1_Line2D), intent(in)    :: self
    real, dimension(:),        intent(in)    :: local_coord
    real, dimension(:,:),      intent(in)    :: elem_coords
    real, dimension(:,:),      intent(inout) :: jacobian
    write(0,*) "LagrangeP1_Line2D::jacobian()"
  end subroutine LagrangeP1_Line2D__jacobian

  subroutine LagrangeP1_Line2D__jacobian_determinant( self, local_coord, elem_coords, jacobian_determinant )
    class(LagrangeP1_Line2D), intent(in)   :: self
    real, dimension(:),       intent(in)   :: local_coord
    real, dimension(:,:),     intent(in)   :: elem_coords
    real, intent(out) :: jacobian_determinant
    jacobian_determinant = 0.
    write(0,*) "LagrangeP1_Line2D::jacobian_determinant()"
  end subroutine LagrangeP1_Line2D__jacobian_determinant


  subroutine LagrangeP1_Quad__init(self)
    class(LagrangeP1_Quad), intent(inout) :: self
    ! Call parent class init
    call self%ShapeFunction_class%init()
    self%shape = 1 ! Quad
    self%dimensionality = 2
    self%nb_nodes = 4
    self%nb_sides = 4
    self%order = 1
    ! Deallocation will happen in base class ShapeFunction__destruct
    allocate( self%local_coords(self%nb_nodes,self%dimensionality) )
    self%local_coords(1,:) = [-1,-1]
    self%local_coords(2,:) = [+1,-1]
    self%local_coords(3,:) = [+1,+1]
    self%local_coords(4,:) = [-1,+1]

    ! print something
    write(0,*) "LagrangeP1_Quad::init"
  end subroutine LagrangeP1_Quad__init  

  subroutine LagrangeP1_Quad__values( self, local_coord, values )
    class(LagrangeP1_Quad), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:), intent(inout) :: values
    values(1) = 0.25 * (1.-local_coord(1)) * (1.-local_coord(2))
    values(2) = 0.25 * (1.+local_coord(1)) * (1.-local_coord(2))
    values(3) = 0.25 * (1.+local_coord(1)) * (1.+local_coord(2))
    values(4) = 0.25 * (1.-local_coord(1)) * (1.+local_coord(2))
    write(0,*) "LagrangeP1_Quad::values"
  end subroutine LagrangeP1_Quad__values

  subroutine LagrangeP1_Quad__grad_values( self, local_coord, grad_values )
    class(LagrangeP1_Quad), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:,:), intent(inout) :: grad_values

    grad_values(1,1) = 0.25 * (-1.+local_coord(2))
    grad_values(1,2) = 0.25 * ( 1.-local_coord(2))
    grad_values(1,3) = 0.25 * ( 1.+local_coord(2))
    grad_values(1,4) = 0.25 * (-1.-local_coord(2))

    grad_values(2,1) = 0.25 * (-1.+local_coord(1))
    grad_values(2,2) = 0.25 * (-1.-local_coord(1))
    grad_values(2,3) = 0.25 * ( 1.+local_coord(1))
    grad_values(2,4) = 0.25 * ( 1.-local_coord(1))

    write(0,*) "LagrangeP1_Quad::grad_values"
  end subroutine LagrangeP1_Quad__grad_values


  subroutine LagrangeP1_Quad2D__init(self)
    class(LagrangeP1_Quad2D), intent(inout) :: self

    ! Specify self%sf
    allocate(LagrangeP1_Quad :: self%sf)
    call self%sf%init()

    ! Call parent class init
    call self%Element%init()    

    self%dimension = 2
    
    ! print something
    write(0,*) "LagrangeP1_Quad2D::init"
  end subroutine LagrangeP1_Quad2D__init
    
    subroutine LagrangeP1_Quad2D__jacobian( self, local_coord, elem_coords, jacobian )
    class(LagrangeP1_Quad2D), intent(in)    :: self
    real, dimension(:),        intent(in)    :: local_coord
    real, dimension(:,:),      intent(in)    :: elem_coords
    real, dimension(:,:),      intent(inout) :: jacobian
    real :: ax,bx,cx,dx,ay,by,cy,dy
    !ax = 0.25*( elem_coords(1,1) + elem_coords(2,1) + elem_coords(3,1) + elem_coords(4,1) )
    bx = 0.25*(-elem_coords(1,1) + elem_coords(2,1) + elem_coords(3,1) - elem_coords(4,1) )
    cx = 0.25*(-elem_coords(1,1) - elem_coords(2,1) + elem_coords(3,1) + elem_coords(4,1) )
    dx = 0.25*( elem_coords(1,1) - elem_coords(2,1) + elem_coords(3,1) - elem_coords(4,1) )
    !ay = 0.25*( elem_coords(1,2) + elem_coords(2,2) + elem_coords(3,2) + elem_coords(4,2) )
    by = 0.25*(-elem_coords(1,2) + elem_coords(2,2) + elem_coords(3,2) - elem_coords(4,2) )
    cy = 0.25*(-elem_coords(1,2) - elem_coords(2,2) + elem_coords(3,2) + elem_coords(4,2) )
    dy = 0.25*( elem_coords(1,2) - elem_coords(2,2) + elem_coords(3,2) - elem_coords(4,2) )
    jacobian(1,1) = bx + dx*local_coord(2)
    jacobian(1,2) = by + dy*local_coord(2)
    jacobian(2,1) = cx + dx*local_coord(1)
    jacobian(2,2) = cy + dy*local_coord(1)
    write(0,*) "LagrangeP1_Quad2D::jacobian()"
  end subroutine LagrangeP1_Quad2D__jacobian

  subroutine LagrangeP1_Quad2D__jacobian_inverse( self, local_coord, elem_coords, jacobian_inverse )
    class(LagrangeP1_Quad2D), intent(in)    :: self
    real, dimension(:),        intent(in)    :: local_coord
    real, dimension(:,:),      intent(in)    :: elem_coords
    real, dimension(:,:),      intent(inout) :: jacobian_inverse
    real, dimension(2,2) :: J
    real :: one_over_Jdet
    call self%jacobian( local_coord, elem_coords, J )
    one_over_Jdet = 1./(J(1,1)*J(2,2)-J(1,2)*J(2,1)) 
    jacobian_inverse(1,1) =  one_over_Jdet*J(2,2)
    jacobian_inverse(1,2) = -one_over_Jdet*J(1,2)
    jacobian_inverse(2,1) = -one_over_Jdet*J(2,1)
    jacobian_inverse(2,2) =  one_over_Jdet*J(1,1)
    write(0,*) "LagrangeP1_Quad2D::jacobian_inverse()"
  end subroutine LagrangeP1_Quad2D__jacobian_inverse
    
  subroutine LagrangeP1_Triag__init(self)
    class(LagrangeP1_Triag), intent(inout) :: self
    ! Call parent class init
    call self%ShapeFunction_class%init()

    self%shape = 2 ! Triag
    self%dimensionality = 2
    self%nb_nodes = 3
    self%nb_sides = 3
    self%order = 1
    ! Deallocation will happen in base class ShapeFunction__destruct
    allocate( self%local_coords(self%nb_nodes,self%dimensionality) )
    self%local_coords(1,:) = [0,0]
    self%local_coords(2,:) = [1,0]
    self%local_coords(3,:) = [0,1]

    ! print something
    write(0,*) "LagrangeP1_Triag::init"
  end subroutine LagrangeP1_Triag__init  

  subroutine LagrangeP1_Triag__values( self, local_coord, values )
    class(LagrangeP1_Triag), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:), intent(inout) :: values
    values(1) = 1. - local_coord(1) - local_coord(2)
    values(2) = local_coord(1)
    values(3) = local_coord(2)
    write(0,*) "LagrangeP1_Triag::values"
  end subroutine LagrangeP1_Triag__values

  subroutine LagrangeP1_Triag__grad_values( self, local_coord, grad_values )
    class(LagrangeP1_Triag), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:,:), intent(inout) :: grad_values
    ! Notice that this is in fact independent of the local_coord
    grad_values(1,1) = -1.
    grad_values(1,2) =  1.
    grad_values(1,3) =  0.

    grad_values(2,1) = -1.
    grad_values(2,2) =  0.
    grad_values(2,3) =  1.

    write(0,*) "LagrangeP1_Triag::grad_values"
  end subroutine LagrangeP1_Triag__grad_values


  subroutine LagrangeP1_Triag2D__init(self)
    class(LagrangeP1_Triag2D), intent(inout) :: self

    ! Specify self%sf
    allocate(LagrangeP1_Triag :: self%sf)
    call self%sf%init()

    ! Call parent class init
    call self%Element%init()    

    self%dimension = 2
    
    ! print something
    write(0,*) "LagrangeP1_Triag2D::init"
  end subroutine LagrangeP1_Triag2D__init

  subroutine LagrangeP1_Triag2D__jacobian( self, local_coord, elem_coords, jacobian )
    class(LagrangeP1_Triag2D), intent(in)    :: self
    real, dimension(:),        intent(in)    :: local_coord
    real, dimension(:,:),      intent(in)    :: elem_coords
    real, dimension(:,:),      intent(inout) :: jacobian
    ! Notice that this is in fact independent of the local_coord
    jacobian(1,1) = elem_coords(2,1) - elem_coords(1,1)
    jacobian(1,2) = elem_coords(2,2) - elem_coords(1,2)
    jacobian(2,1) = elem_coords(3,1) - elem_coords(1,1)
    jacobian(2,2) = elem_coords(3,2) - elem_coords(1,2)
    write(0,*) "LagrangeP1_Triag2D::jacobian()"
  end subroutine LagrangeP1_Triag2D__jacobian

  subroutine LagrangeP1_Triag2D__jacobian_inverse( self, local_coord, elem_coords, jacobian_inverse )
    class(LagrangeP1_Triag2D), intent(in)    :: self
    real, dimension(:),        intent(in)    :: local_coord
    real, dimension(:,:),      intent(in)    :: elem_coords
    real, dimension(:,:),      intent(inout) :: jacobian_inverse
    real, dimension(2,2) :: J
    real :: one_over_Jdet
    call self%jacobian( local_coord, elem_coords, J )
    one_over_Jdet = 1./(J(1,1)*J(2,2)-J(1,2)*J(2,1)) 
    jacobian_inverse(1,1) =  one_over_Jdet*J(2,2)
    jacobian_inverse(1,2) = -one_over_Jdet*J(1,2)
    jacobian_inverse(2,1) = -one_over_Jdet*J(2,1)
    jacobian_inverse(2,2) =  one_over_Jdet*J(1,1)
    write(0,*) "LagrangeP1_Triag2D::jacobian_inverse()"
  end subroutine LagrangeP1_Triag2D__jacobian_inverse
 
end module lagrangep1_module
