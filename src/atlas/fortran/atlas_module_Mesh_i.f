! (C) Copyright 2013-2014 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(object_type) :: Mesh_type

! Purpose :
! -------
!   *Mesh* : Container type holding an entire mesh

! Methods :
! -------
!   add_function_space : Add a new function space
!   function_space : Access the function space with given name

! Author :
! ------
!   20-Nov-2013 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure :: create_function_space => Mesh__create_function_space
  procedure :: function_space => Mesh__function_space
END TYPE Mesh_type

interface new_Mesh
  module procedure new_Mesh
end interface
!------------------------------------------------------------------------------

