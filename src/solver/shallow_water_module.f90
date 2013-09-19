
! module describing a shallow_water model
module shallow_water_module
  use grid_module
  use model_module
  use mpdata_module
  implicit none

  type, public, extends(State_class) :: ShallowWaterState
    class(Field_class),  pointer :: D, Q, H0
  contains
    procedure, pass :: init => ShallowWaterState__init
  end type ShallowWaterState


  type, public, extends(Model_class) :: ShallowWaterModel
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
    class(Field_class), private, pointer :: grad_D ! gradient of depth
    class(Field_class), private, pointer :: R, Q0, D0, Q, D, V

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

 function new_ShallowWaterState(name, function_space) result(state)
    character(len=*), intent(in)              :: name
    class(FunctionSpace_class), pointer, intent(in) :: function_space
    class(State_class), pointer :: state
    allocate( ShallowWaterState :: state )
    call state%init(name,function_space)
  end function new_ShallowWaterState

  subroutine ShallowWaterState__init(self,name,function_space)
    class(ShallowWaterState), intent(inout) :: self
    character(len=*), intent(in) :: name
    class(FunctionSpace_class), intent(in), target :: function_space
    call self%State_class%init(name,function_space)

    write(0,*) "ShallowWaterState::init(",name,",",function_space%name,")"  
    !self%D  => self%function_space%add_scalar_field("depth")
    !self%Q  => self%function_space%add_vector_field("momentum")
    !self%H0 => self%function_space%add_scalar_field("orography")
    self%D  => new_ScalarField("depth",function_space)
    self%Q  => new_VectorField("momentum",function_space)
    self%H0 => new_ScalarField("orography",function_space)
    call self%add_field(self%D)
    call self%add_field(self%Q)
    call self%add_field(self%H0)
  end subroutine ShallowWaterState__init


! ------------------------------------------------------------------------------------
!                              ShallowWaterModel subroutines
!-------------------------------------------------------------------------------------

  subroutine ShallowWaterModel__init(self,grid)
    class(ShallowWaterModel), intent(inout), target :: self
    class(Grid_class), intent(in), target :: grid
    class(FunctionSpace_class), pointer :: vertices
    class(Field_class), pointer :: dual_volume
    class(Field_class), pointer :: dual_face_normal

    call self%Model_class%init(grid)
    write(0,*) "ShallowWaterModel::init(grid)"
 
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
    class(FunctionSpace_class), pointer :: vertices
    class(Field_class), pointer :: dual_volume
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
    class(Field_class) , pointer :: D, Q
    real, dimension(:,:), pointer :: nodes
    integer :: inode
    integer :: nb_nodes

    integer :: ir,ip
    real    :: aaa0,zk,om,ph0,g,ath,bth,cth,pi,x,y,th,cor

    ! statement-functions, before first statement
    ATH(TH) = om*0.5*(self%f0+om)*(cos(TH))**2 &
      & +0.25*zk**2*(cos(TH))**(2*ir)*( (ir+1)*(cos(TH))**2 &
      & +real(2*ir**2-ir-2)-2.*ir**2/(cos(TH))**2 )
    BTH(TH) = (self%f0+2.*om)*zk/real((ir+1)*(ir+2))*(cos(TH))**ir &
      & *( real(ir**2+2*ir+2)-((ir+1)*cos(TH))**2 )
    CTH(TH) = 0.25*zk**2*(cos(TH))**(2*ir)*( real(ir+1)*(cos(TH))**2 &
      & -real(ir+2) )  

    om   = 7.848E-6
    zk   = 7.848E-6
    ir   = 4
    ph0  = 78.4E3
    aaa0 = 0.

    D => self%state%field("depth")
    Q => self%state%field("momentum")

    nodes    => self%state%function_space%grid%nodes
    nb_nodes = self%state%function_space%grid%nb_nodes

      do inode=1,nb_nodes
        x=nodes(inode,1)
        y=nodes(inode,2)
        if(x == 2.*self%pi) x=0.
          cor=self%f0*sin(y)
          Q%array(inode,1) =  self%radius*OM*cos(y)+self%radius*ZK*cos(IR*x) *(cos(y))**(IR-1)*(IR*(sin(y))**2-(cos(y))**2)
          Q%array(inode,2) = -self%radius*ZK*IR*(cos(y))**(IR-1)*sin(y)*sin(IR*x)
          D%array(inode,1) = (ph0+self%radius**2*ATH(y)+self%radius**2*BTH(y)*cos(IR*x)+self%radius**2*CTH(y)*cos(2.*IR*x)) &
              & /(self%grav)
          D%array(inode,1) = max(aaa0,D%array(inode,1))
          Q%array(inode,1) = Q%array(inode,1) * D%array(inode,1)
          Q%array(inode,2) = Q%array(inode,2) * D%array(inode,1)
          if(y == 0.5*self%pi) Q%array(inode,1)=0.
          if(y ==-0.5*self%pi) Q%array(inode,1)=0.
      end do

  end subroutine ShallowWaterModel__set_state_rossby_haurwitz
  

! ==========================================================================
! MPDATA Solver subroutines
! ==========================================================================

  function new_ShallowWaterSolver(model) result(solver)
    class(Model_class), pointer :: model
    class(Solver_class), pointer :: solver
    allocate( ShallowWaterSolver :: solver )
    write(0,*) "new_ShallowWaterSolver"
    call solver%init(model)
  end function new_ShallowWaterSolver

  subroutine ShallowWaterSolver__init(self,model)
    class(ShallowWaterSolver), intent(inout) :: self
    class(Model_class), intent(in), target :: model
    class(FunctionSpace_class), pointer :: vertices
    class(FunctionSpace_class), pointer :: faces

    call self%MPDATA_Solver%init(model)
    self%state => self%model%state
    write(0,*) "ShallowWaterSolver::init(model)"  
    vertices    => self%model%grid%function_space("vertices")
    self%D      => self%state%field("depth")
    self%Q      => self%state%field("momentum")

    self%grad_D => new_VectorField("depth_gradient",vertices)
    self%R      => new_VectorField("rhs_forcing",vertices)
    self%D0     => new_ScalarField("depth_prev_iter",vertices)
    self%Q0     => new_VectorField("momentum_prev_iter",vertices)
    self%V      => new_VectorField("advective_velocity",vertices)

    call self%state%add_field(self%R)
    call self%state%add_field(self%grad_D)

  end subroutine ShallowWaterSolver__init

  subroutine ShallowWaterSolver__backup_solution(self)
    class(ShallowWaterSolver), intent(inout) :: self
    self%Q0%array = self%Q%array
    self%D0%array = self%D%array
  end subroutine ShallowWaterSolver__backup_solution

  subroutine ShallowWaterSolver__add_forcing_to_solution(self)
    class(ShallowWaterSolver), intent(inout) :: self
    integer :: inode
    do inode=1,self%model%grid%nb_nodes
      self%Q%array(inode,1) = self%Q%array(inode,1) + 0.5*self%dt*self%R%array(inode,1)
      self%Q%array(inode,2) = self%Q%array(inode,2) + 0.5*self%dt*self%R%array(inode,2)
    end do
  end subroutine ShallowWaterSolver__add_forcing_to_solution

  subroutine ShallowWaterSolver__compute_advective_velocities(self)
    class(ShallowWaterSolver), intent(inout)    :: self
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
      !self%V%array(inode,1) = Qx/D * r
      !self%V%array(inode,2) = Qy/D * r*cos(y)
      !self%V%array(inode,1) = 50.*r*cos(y)
      !self%V%array(inode,2) = 0.     
    end do
  end subroutine ShallowWaterSolver__compute_advective_velocities

  subroutine ShallowWaterSolver__advect_solution(self)
    class(ShallowWaterSolver), intent(inout)    :: self
    integer :: inode, p1, p2
    call self%mpdata_gauge( self%D%array(:,1), self%V%array, .False. )
    call self%mpdata_gauge( self%Q%array(:,1), self%V%array, .True. )
    call self%mpdata_gauge( self%Q%array(:,2), self%V%array, .True. )

    ! remove noise from periodic boundary.. why is it even there in the mesh?
    do inode=1,self%model%grid%nb_periodic_nodes
      p1 = self%model%grid%periodic_nodes(inode,1)
      p2 = self%model%grid%periodic_nodes(inode,2)
      self%D%array(p2,1) = self%D%array(p1,1)
      self%Q%array(p2,1) = self%Q%array(p1,1)
      self%Q%array(p2,2) = self%Q%array(p1,2)
    end do
  end subroutine ShallowWaterSolver__advect_solution

  subroutine ShallowWaterSolver__compute_Rn(self)
    class(ShallowWaterSolver), intent(inout)    :: self
    real :: f0, f, x, y, r, g, Qx, Qy, D, dDdx, dDdy, sin_y, cos_y, eps, vol, hx, hy
    integer :: inode
    
    eps = 1e-6
    r   = 6371.22e+03
    f0  = 1.4584e-04 !coriolis parameter (=2xearth's omega)
    g   = 9.80616


    call self%compute_gradient( self%D%array(:,1), self%grad_D%array, .False. )
    
    do inode=1,self%model%grid%nb_nodes
      y     = self%model%grid%nodes(inode,2)
      sin_y = sin(y)
      cos_y = max( eps, cos(y) )
      Qx    = self%Q%array(inode,1)
      Qy    = self%Q%array(inode,2)
      D     = max( eps, self%D%array(inode,1) )
      dDdx  = self%grad_D%array(inode,1)
      dDdy  = self%grad_D%array(inode,2)
      vol   = self%vol%array(inode,1)
      hx    = r*cos_y
      hy    = r
      f     = f0 * sin_y
      self%R%array(inode,1)   = -g*D*dDdx*hy/vol !+ f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
      self%R%array(inode,2)   = -g*D*dDdy*hx/vol !- f*Qx - sin_y/(r*cos_y)*Qx*Qx/D
    end do
  end subroutine ShallowWaterSolver__compute_Rn


  subroutine ShallowWaterSolver__implicit_solve(self)
    class(ShallowWaterSolver), intent(inout)    :: self
    real :: f0, f, y, r, g, Rx, Ry, Qx, Qy, D, dDdx, dDdy, hx, hy, sin_y, cos_y, eps, Qx_adv, Qy_adv, Rx_exp, Ry_exp, vol
    integer :: inode, m

    eps = 1e-10
    r   = 6371.22e+03
    f0  = 1.4584e-04 !coriolis parameter (=2xearth's omega)
    g   = 9.80616

    ! D is already up to date at time level (n+1), just by MPDATA advection
    call self%compute_gradient( self%D%array(:,1), self%grad_D%array, .False. )
    
    do inode=1,self%model%grid%nb_nodes
      y     = self%model%grid%nodes(inode,2)
      sin_y = sin(y)
      cos_y = max( eps, cos(y) )
      D     = max( eps, self%D%array(inode,1) )
      dDdx  = self%grad_D%array(inode,1)
      dDdy  = self%grad_D%array(inode,2)
      vol   = self%vol%array(inode,1)
      hx    = r*cos_y
      hy    = r
      f     = f0 * sin_y

      Rx_exp = -g*D*dDdx*hy/vol
      Ry_exp = -g*D*dDdy*hx/vol

      Qx = self%Q%array(inode,1)
      Qy = self%Q%array(inode,2)
      Qx_adv = Qx
      Qy_adv = Qy

      do m=1,3 ! Three iterations at most is enough to converge
        Rx = Rx_exp + f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
        Ry = Ry_exp - f*Qx - sin_y/(r*cos_y)*Qx*Qx/D
      
        Qx = Qx_adv + 0.5*self%dt*Rx
        Qy = Qy_adv + 0.5*self%dt*Ry
      end do

      self%Q%array(inode,1) = Qx
      self%Q%array(inode,2) = Qy

      Rx = Rx_exp + f*Qy + sin_y/(r*cos_y)*Qx*Qy/D
      Ry = Ry_exp - f*Qx - sin_y/(r*cos_y)*Qx*Qx/D      
      self%R%array(inode,1) = Rx
      self%R%array(inode,2) = Ry
    end do

  end subroutine ShallowWaterSolver__implicit_solve

end module shallow_water_module
