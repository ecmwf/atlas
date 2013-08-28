
! module describing a shallow_water model
module shallow_water_module
  use grid_module
  use state_module
  implicit none

  type, public, extends(State) :: ShallowWaterState
    type(Field),  pointer :: D, Q
  contains
    procedure, pass :: init => ShallowWaterState__init
  end type ShallowWaterState

contains

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
    self%D => self%function_space%add_scalar_field("depth")
    self%Q => self%function_space%add_vector_field("momentum")
    call self%add_field(self%D)
    call self%add_field(self%Q)

  end subroutine ShallowWaterState__init


  subroutine init_state_rossby_haurwitz(S)
    class(State), intent(in) :: S
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
    g    = 9.80616
    aaa0 = 0.
    a    = 6371.22e+03 ! Earth radius
    pi   = acos(-1.)
    f0   = 1.4584e-04 !coriolis parameter (=2xearth's omega)

    select type (S)
    type is (ShallowWaterState)
      D => S%D
      Q => S%Q

      nodes => S%function_space%g%nodes
      nb_nodes = S%function_space%g%nb_nodes

      do inode=1,nb_nodes
        x=nodes(inode,1)
        y=nodes(inode,2)
        if(x == 2.*pi) x=0.
          cor=f0*sin(y)
          Q%array(inode,1) =  a*OM*cos(y)+a*ZK*cos(IR*x) *(cos(y))**(IR-1)*(IR*(sin(y))**2-(cos(y))**2)
          Q%array(inode,2) = -a*ZK*IR*(cos(y))**(IR-1)*sin(Y)*sin(IR*x)
          D%array(inode,1) = ( PH0+a**2*ATH(y)+a**2*BTH(y)*cos(IR*x)+a**2*CTH(y)*cos(2.*IR*x) ) / g
          D%array(inode,1) = max(aaa0,D%array(inode,1))
          Q%array(inode,1) = Q%array(inode,1) * D%array(inode,1)
          Q%array(inode,2) = Q%array(inode,2) * D%array(inode,1)
          if(y == 0.5*pi) Q%array(inode,1)=0.
          if(y ==-0.5*pi) Q%array(inode,1)=0.
      enddo
    end select

  end subroutine init_state_rossby_haurwitz

end module shallow_water_module
