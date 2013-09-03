module elements_module
  implicit none
  private
  
  ! Shapefunction class
  ! -------------------
  ! Provides local-coordinate (XI,ETA,ZETA) functionality
  ! - interpolation  Q(XI,ETA,ZETA) = sum( L_i(XI,ETA,ZETA) * Q_i )
  ! - gradient       [dQ/dXI, dQ/dETA, dQ/dZETA]
  ! - divergence     dQ/dXI + dQ/dETA + dQ/dZETA
  type, public :: ShapeFunction_class
    integer :: shape ! 1=Quad, 2=Triag
    integer :: dimensionality
    integer :: nb_nodes
    integer :: nb_sides
    integer :: order
    real, dimension(:,:), allocatable :: local_coords
  contains
    procedure, public :: init          => ShapeFunction__init
    procedure, public :: destruct      => ShapeFunction__destruct
    procedure, public :: values        => ShapeFunction__values
    procedure, public :: grad_values   => ShapeFunction__grad_values

  end type ShapeFunction_class

  ! Element class
  ! -------------------
  ! Provides link between shapefunction and real coordinate system
  ! A Quad element could exist in the 2D world or in a 3D world as on a sphere,
  ! While the shapefunction Quad will be the same in its local coordinate system
  ! - jacobian transformation matrix
  type, public :: Element
    class(ShapeFunction_class), allocatable :: sf
    integer :: dimension
    integer :: dimensionality
    integer :: nb_nodes
    integer :: nb_sides
  contains
    procedure, public :: init     => Element__init
    procedure, public :: destruct => Element__destruct
    procedure, public :: jacobian => Element__jacobian
    procedure, public :: jacobian_inverse => Element__jacobian_inverse

    procedure, public :: print    => Element__print
  end type Element

contains
  
  ! ShapeFunction::init
  ! -------------------
  ! Constructor for ShapeFunction base class
  subroutine ShapeFunction__init(self)
    class(ShapeFunction_class), intent(inout) :: self
    write(0,*) "ShapeFunction::init"
  end subroutine ShapeFunction__init
  
  ! ShapeFunction::destruct
  ! -------------------
  ! Destructor for ShapeFunction base class
  subroutine ShapeFunction__destruct(self)
    class(ShapeFunction_class) :: self
    write(0,*) "ShapeFunction::destruct"
    if (allocated(self%local_coords)) deallocate(self%local_coords)
  end subroutine ShapeFunction__destruct

  ! ShapeFunction::values
  ! ---------------------
  ! Compute the interpolation values for a given local coordinate
  subroutine ShapeFunction__values( self, local_coord, values )
    class(ShapeFunction_class), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:), intent(inout) :: values

    write(0,*) "ShapeFunction::values"
  end subroutine ShapeFunction__values

  ! ShapeFunction::grad_values
  ! --------------------------
  ! Compute the interpolation values for a given local coordinate
  subroutine ShapeFunction__grad_values( self, local_coord, grad_values )
    class(ShapeFunction_class), intent(in) :: self
    real, dimension(:), intent(in) :: local_coord
    real, dimension(:,:), intent(inout) :: grad_values

    write(0,*) "ShapeFunction::grad_values"
  end subroutine ShapeFunction__grad_values

  ! Element::init
  ! -------------------
  ! Constructor for Element base class
  subroutine Element__init(self)
    class(Element), intent(inout) :: self
    self%dimensionality = self%sf%dimensionality
    self%nb_nodes = self%sf%nb_nodes
    self%nb_sides = self%sf%nb_sides  
    write(0,*) "Element::init"
  end subroutine Element__init
  
  ! Element::destruct
  ! -------------------
  ! Destructor for Element base class
  subroutine Element__destruct(self)
    class(Element), intent(inout) :: self
    write(0,*) "Element::destruct"
    call self%sf%destruct()
    if( allocated(self%sf) ) deallocate(self%sf)
  end subroutine Element__destruct
  
  ! Element::print
  ! --------------
  ! Print all information from the element type
  subroutine Element__print(self)
    class(Element), intent(inout) :: self
    integer :: n
    select case(self%sf%shape)
      case(1)
    write(0,*) "shape           = QUAD"
      case(2)
    write(0,*) "shape           = TRIAG"
    end select
    write(0,*) "order           = ",self%sf%order
    write(0,*) "dimensionality  = ",self%dimensionality
    write(0,*) "dimension       = ",self%dimension
    write(0,*) "nb_nodes        = ",self%nb_nodes
    write(0,*) "nb_sides        = ",self%nb_sides
    write(0,*) "local_coords    = "
    do n = 1,self%nb_nodes
      write(0,*) "  ",self%sf%local_coords(n,:)
    end do
  end subroutine Element__print

  ! Element::jacobian
  ! -----------------
  ! Compute the jacobian
  ! [ dx/dxi   dy/dxi   dz/dxi   ]
  ! [ dx/deta  dy/deta  dz/deta  ]
  ! [ dx/dzeta dy/dzeta dz/dzeta ]
  subroutine Element__jacobian( self, local_coord, elem_coords, jacobian )
    class(Element), intent(in)          :: self
    real, dimension(:), intent(in)      :: local_coord
    real, dimension(:,:), intent(in)    :: elem_coords
    real, dimension(:,:), intent(inout) :: jacobian
    write(0,*) "Element::jacobian()"
  end subroutine Element__jacobian

   subroutine Element__jacobian_inverse( self, local_coord, elem_coords, jacobian_inverse )
    class(Element), intent(in)    :: self
    real, dimension(:),        intent(in)    :: local_coord
    real, dimension(:,:),      intent(in)    :: elem_coords
    real, dimension(:,:),      intent(inout) :: jacobian_inverse
    write(0,*) "Element::jacobian_inverse()"
  end subroutine Element__jacobian_inverse
 
end module elements_module