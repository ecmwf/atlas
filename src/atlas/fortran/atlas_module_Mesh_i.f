! (C) Copyright 2013-2014 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_Mesh

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
  procedure, private :: create_function_space_nodes => Mesh__create_function_space_nodes
  procedure, private :: create_function_space_shape => Mesh__create_function_space_shape
  generic :: create_function_space => create_function_space_nodes, create_function_space_shape
  procedure :: function_space => Mesh__function_space
END TYPE atlas_Mesh

interface new_atlas_Mesh
  module procedure new_atlas_Mesh
end interface
!------------------------------------------------------------------------------

