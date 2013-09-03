
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

  type, public, extends(MPDATA_Solver) :: ShallowWaterSolver
    type(Field), private, pointer :: grad_D ! gradient of depth
    type(Field), private, pointer :: R, Rex, Qadv, Q0, D0, Q, D, V

  contains
    procedure, public,  pass :: init => ShallowWaterSolver__init
    procedure, pass :: compute_Rn   => ShallowWaterSolver__compute_Rn
    procedure, pass :: implicit_solve => ShallowWaterSolver__implicit_solve
    procedure, pass :: compute_advective_velocities => ShallowWaterSolver__compute_advective_velocities
    procedure, pass :: backup_solution => ShallowWaterSolver__backup_solution
    procedure, pass :: add_forcing_to_solution => ShallowWaterSolver__add_forcing_to_solution
    procedure, pass :: advect_solution => ShallowWaterSolver__advect_solution
    
  end type ShallowWaterSolver
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
    class(Grid_class), intent(in), target :: g
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
    self%solver => new_ShallowWaterSolver(self%ptr)
    self%solver%state => self%state
  end subroutine ShallowWaterModel__init

  subroutine ShallowWaterModel__compute_metrics(self)
    class(ShallowWaterModel), intent(inout) :: self
    class(FunctionSpace), pointer :: vertices
    type(Field), pointer :: dual_volume
    real :: y, hx, hy, jac
    integer :: inode

    vertices => self%grid%function_space("vertices")
    dual_volume => vertices%field("dual_volume")

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
  

! ==========================================================================
! MPDATA Solver subroutines
! ==========================================================================

  function new_ShallowWaterSolver(model_) result(solver_)
    class(Model), pointer :: model_
    class(Solver), pointer :: solver_
    allocate( ShallowWaterSolver :: solver_ )
    write(0,*) "new_ShallowWaterSolver"
    call solver_%init(model_)
  end function new_ShallowWaterSolver

  subroutine ShallowWaterSolver__init(self,model_)
    class(ShallowWaterSolver), intent(inout) :: self
    class(Model), intent(in), target :: model_
    class(FunctionSpace), pointer :: vertices
    class(FunctionSpace), pointer :: faces

    call self%MPDATA_Solver%init(model_)
    
    
    write(0,*) "ShallowWaterSolver::init(model)"  
    vertices => self%model%grid%function_space("vertices")
    self%D => vertices%field("depth")
    self%Q => vertices%field("momentum")
    self%grad_D => vertices%add_vector_field("depth_gradient")
    self%R => vertices%add_vector_field("rhs_forcing")
    self%Rex => vertices%add_vector_field("rhs_forcing_exact_part")
    self%Qadv => vertices%add_vector_field("momentum_advected")
    self%D0 => vertices%add_vector_field("depth_prev_iter")
    self%Q0 => vertices%add_vector_field("momentum_prev_iter")
    self%V => vertices%add_vector_field("advective_velocity")

  end subroutine ShallowWaterSolver__init

  subroutine ShallowWaterSolver__backup_solution(self)
    class(ShallowWaterSolver), intent(in) :: self
    self%Q0%array = self%Q%array
    self%D0%array = self%D%array
  end subroutine ShallowWaterSolver__backup_solution

  subroutine ShallowWaterSolver__add_forcing_to_solution(self)
    class(ShallowWaterSolver), intent(in) :: self
    integer :: inode
    do inode=1,self%model%grid%nb_nodes
      self%Q%array(inode,1) = self%Q%array(inode,1) + 0.5*self%dt*self%R%array(inode,1)
      self%Q%array(inode,2) = self%Q%array(inode,2) + 0.5*self%dt*self%R%array(inode,2)
    end do
  end subroutine ShallowWaterSolver__add_forcing_to_solution

  subroutine ShallowWaterSolver__compute_advective_velocities(self)
    class(ShallowWaterSolver), intent(in)    :: self
    real :: y, r, Qx, Qy, Q0x, Q0y, D, D0, eps
    integer :: inode
    
    eps = 1e-10
    r = 6371.22e+03
   
    do inode=1,self%model%grid%nb_nodes
      y     = self%model%grid%nodes(inode,2)
      Qx    = self%Q%array(inode,1)
      Qy    = self%Q%array(inode,2)
      D     = max( eps, self%D%array(inode,1) )
      Q0x   = self%Q0%array(inode,1)
      Q0y   = self%Q0%array(inode,2)
      D0    = max( eps, self%D0%array(inode,1) )
      
      ! this really computes V = G*contravariant_velocity, 
      ! with    G=hx*hy,
      !         physical_velocity = dotproduct( [hx,hy] , contravariant_velocity )
      ! V = (hx*hy) * [u/hx, v/hy] = [u/hy, v/hx]
      ! and hx = r*cos(y)  ,  hy = r
      ! and Q = [ D*u , D*v ]
      self%V%array(inode,1) = ( 1.5*Qx/D - 0.5*Q0x/D0 )*r
      self%V%array(inode,2) = ( 1.5*Qy/D - 0.5*Q0y/D0 )*r*cos(y)
      
    end do
  end subroutine ShallowWaterSolver__compute_advective_velocities

  subroutine ShallowWaterSolver__advect_solution(self)
    class(ShallowWaterSolver), intent(in)    :: self
    call self%mpdata_gage( self%D%array(:,1), self%V%array, self%D%array(:,1)    )
    call self%mpdata_gage( self%Q%array(:,1), self%V%array, self%Qadv%array(:,1) )
    call self%mpdata_gage( self%Q%array(:,2), self%V%array, self%Qadv%array(:,2) )
  end subroutine ShallowWaterSolver__advect_solution

  subroutine ShallowWaterSolver__compute_Rn(self)
    class(ShallowWaterSolver), intent(in)    :: self
    real :: f0, f, x, y, r, g, Qx, Qy, D, dDdx, sin_y, cos_y, eps
    integer :: inode
    
    eps = 1e-10
    r   = 6371.22e+03
    f0  = 1.4584e-04 !coriolis parameter (=2xearth's omega)
    g   = 9.80616


    call self%compute_gradient( self%D%array(:,1), self%grad_D%array )
    
    do inode=1,self%model%grid%nb_nodes
      y     = self%model%grid%nodes(inode,2)
      sin_y = sin(y)
      cos_y = max( eps, cos(y) )
      Qx    = self%Q%array(inode,1)
      Qy    = self%Q%array(inode,2)
      D     = max( eps, self%D%array(inode,1) )
      dDdx  = self%grad_D%array(inode,1)
      f     = f0 * sin_y
      
      self%Rex%array(inode,1) = -g/(r*cos_y)*D*dDdx
      self%Rex%array(inode,2) = -g/r*D*dDdx
      self%R%array(inode,1)   = self%Rex%array(inode,1) + f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
      self%R%array(inode,2)   = self%Rex%array(inode,2) - f*Qx - sin_y/(r*cos_y)*Qx*Qx/D
    end do
  end subroutine ShallowWaterSolver__compute_Rn


  subroutine ShallowWaterSolver__implicit_solve(self)
    class(ShallowWaterSolver), intent(in)    :: self
    real :: f0, f, y, r, g, Qx, Qy, D, dDdx, sin_y, cos_y, eps
    integer :: inode, m

    eps = 1e-10
    r   = 6371.22e+03
    f0  = 1.4584e-04 !coriolis parameter (=2xearth's omega)
    g   = 9.80616

    ! D is already up to date at time level (n+1), just by MPDATA advection
    call self%compute_gradient( self%D%array(:,1), self%grad_D%array )
    
    do inode=1,self%model%grid%nb_nodes
      y     = self%model%grid%nodes(inode,2)
      Qx    = self%Q%array(inode,1)
      Qy    = self%Q%array(inode,2)
      cos_y = max( eps, cos(y) )
      D     = self%D%array(inode,1)
      dDdx  = self%grad_D%array(inode,1)
      self%Rex%array(inode,1) = -g/(r*cos_y)*D*dDdx
      self%Rex%array(inode,2) = -g/r*D*dDdx
    end do

    do inode=1,self%model%grid%nb_nodes
      y     = self%model%grid%nodes(inode,2)
      sin_y = sin(y)
      cos_y = max( eps, cos(y) )
      D     = max( eps, self%D%array(inode,1) )
      f     = f0 * sin_y

      do m=1,3 ! Three iterations at most is enough to converge
        Qx = self%Q%array(inode,1)
        Qy = self%Q%array(inode,2)

        self%R%array(inode,1) = self%Rex%array(inode,1) + f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
        self%R%array(inode,2) = self%Rex%array(inode,2) - f*Qx - sin_y/(r*cos_y)*Qx*Qx/D
      
        self%Q%array(inode,1) = self%Qadv%array(inode,1) + 0.5*self%dt*self%R%array(inode,1)
        self%Q%array(inode,2) = self%Qadv%array(inode,2) + 0.5*self%dt*self%R%array(inode,2)
      end do
    
    end do
  end subroutine ShallowWaterSolver__implicit_solve
  

end module shallow_water_module
