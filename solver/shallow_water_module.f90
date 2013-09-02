
! module describing a shallow_water model
module shallow_water_module
  use grid_module
  use model_module
  use mpdata_module
  implicit none

  type, public, extends(State) :: ShallowWaterState
    type(Field),  pointer :: D, Q, H0
  contains
    procedure, pass :: init => ShallowWaterState__init
  end type ShallowWaterState

  type, public, extends(Model) :: ShallowWaterModel
    real :: radius
    real :: grav
    real :: f0
    real :: pi
  contains
    procedure, pass :: init => ShallowWaterModel__init
    procedure, pass :: compute_metrics => ShallowWaterModel__compute_metrics
    procedure, pass :: set_state_rossby_haurwitz => ShallowWaterModel__set_state_rossby_haurwitz
  end type ShallowWaterModel

contains

! ------------------------------------------------------------------------------------
!                              ShallowWaterState subroutines
!-------------------------------------------------------------------------------------

 function new_ShallowWaterState(name, function_space) result(state_)
    character(len=*), intent(in)              :: name
    class(FunctionSpace), pointer, intent(in) :: function_space
    class(State), pointer :: state_
    allocate( ShallowWaterState :: state_ )
    call state_%init(name,function_space)
  end function new_ShallowWaterState

  subroutine ShallowWaterState__init(self,name,function_space)
    class(ShallowWaterState), intent(inout) :: self
    character(len=*), intent(in) :: name
    class(FunctionSpace), intent(in), target :: function_space
    call self%State%init(name,function_space)

    write(0,*) "ShallowWaterState::init(",name,",",function_space%name,")"  
    self%D  => self%function_space%add_scalar_field("depth")
    self%Q  => self%function_space%add_vector_field("momentum")
    self%H0 => self%function_space%add_scalar_field("orography")
    call self%add_field(self%D)
    call self%add_field(self%Q)
    call self%add_field(self%H0)
  end subroutine ShallowWaterState__init


! ------------------------------------------------------------------------------------
!                              ShallowWaterModel subroutines
!-------------------------------------------------------------------------------------

  subroutine ShallowWaterModel__init(self,g)
    class(ShallowWaterModel), intent(inout), target :: self
    class(Grid), intent(in), target :: g
    class(FunctionSpace), pointer :: vertices
    type(Field), pointer :: dual_volume
    type(Field), pointer :: dual_face_normal

    call self%Model%init(g)
    write(0,*) "ShallowWaterModel::init()"
 
    vertices => self%grid%function_space("vertices")

     ! Coordinate transformation correction for dual_volume on sphere
    self%radius = 6371.22e+03
    self%f0     = 1.4584e-04 !coriolis parameter (=2xearth's omega)
    self%grav   = 9.80616
    self%pi     = acos(-1.)

    ! When changing above constants,
    ! this function should be called AFTER INSTEAD of HERE 
    call self%compute_metrics()

    self%state => new_ShallowWaterState("state", vertices )
    self%solver => new_MPDATA_Solver(self%ptr)
    self%solver%state => self%state
  end subroutine ShallowWaterModel__init

  subroutine ShallowWaterModel__compute_metrics(self)
    class(ShallowWaterModel), intent(inout) :: self
    class(FunctionSpace), pointer :: vertices
    type(Field), pointer :: dual_volume
    real :: y, hx, hy, jac
    integer :: inode

    vertices => self%grid%function_space("vertices")
    dual_volume => vertices%field('dual_volume')

    do inode=1,self%grid%nb_nodes
      y=self%grid%nodes(inode,2)
      hx=self%radius*cos(y)
      hy=self%radius
     jac=hx*hy
     dual_volume%array(inode,1) = dual_volume%array(inode,1)*jac
    enddo
  end subroutine ShallowWaterModel__compute_metrics


  subroutine ShallowWaterModel__set_state_rossby_haurwitz(self)
    class(ShallowWaterModel), intent(inout) :: self
    type(Field),  pointer                   :: D, Q
    real, dimension(:,:), pointer :: nodes
    integer :: inode
    integer :: nb_nodes

    integer :: ir,ip
    real    :: aaa0,zk,om,ph0,g,ath,bth,cth,pi,f0,x,y,th,a,cor

    ! statement-functions, before first statement
    ATH(TH) = om*0.5*(f0+om)*(cos(TH))**2 &
      & +0.25*zk**2*(cos(TH))**(2*ir)*( (ir+1)*(cos(TH))**2 &
      & +real(2*ir**2-ir-2)-2.*ir**2/(cos(TH))**2 )
    BTH(TH) = (f0+2.*om)*zk/real((ir+1)*(ir+2))*(cos(TH))**ir &
      & *( real(ir**2+2*ir+2)-((ir+1)*cos(TH))**2 )
    CTH(TH) = 0.25*zk**2*(cos(TH))**(2*ir)*( real(ir+1)*(cos(TH))**2 &
      & -real(ir+2) )  

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3
    aaa0 = 0.

    select type ( state_ =>self%state )
    type is (ShallowWaterState)
      D => state_%D
      Q => state_%Q

      nodes    => state_%function_space%g%nodes
      nb_nodes = state_%function_space%g%nb_nodes

      do inode=1,nb_nodes
        x=nodes(inode,1)
        y=nodes(inode,2)
        if(x == 2.*self%pi) x=0.
          cor=self%f0*sin(y)
          Q%array(inode,1) =  self%radius*OM*cos(y)+self%radius*ZK*cos(IR*x) *(cos(y))**(IR-1)*(IR*(sin(y))**2-(cos(y))**2)
          Q%array(inode,2) = -self%radius*ZK*IR*(cos(y))**(IR-1)*sin(Y)*sin(IR*x)
          D%array(inode,1) = (PH0+self%radius**2*ATH(y)+self%radius**2*BTH(y)*cos(IR*x)+self%radius**2*CTH(y)*cos(2.*IR*x)) &
              & /(self%grav)
          D%array(inode,1) = max(aaa0,D%array(inode,1))
          Q%array(inode,1) = Q%array(inode,1) * D%array(inode,1)
          Q%array(inode,2) = Q%array(inode,2) * D%array(inode,1)
          if(y == 0.5*self%pi) Q%array(inode,1)=0.
          if(y ==-0.5*self%pi) Q%array(inode,1)=0.
      end do
    end select

  end subroutine ShallowWaterModel__set_state_rossby_haurwitz

end module shallow_water_module
