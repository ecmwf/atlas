! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_Mesh

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
#if ! DEPRECATE_OLD_FUNCTIONSPACE
  procedure, private :: create_function_space_nodes => Mesh__create_function_space_nodes
  procedure, private :: create_function_space_shape => Mesh__create_function_space_shape
  generic :: create_function_space => create_function_space_nodes, create_function_space_shape
  procedure :: function_space => Mesh__function_space
#endif
  procedure :: create_nodes => Mesh__create_nodes
  procedure :: nodes => Mesh__nodes
  procedure :: cells => Mesh__cells
  procedure :: edges => Mesh__edges
  procedure, public :: delete => atlas_Mesh__delete
  procedure, public :: copy => atlas_Mesh__copy
END TYPE atlas_Mesh

interface atlas_Mesh
  module procedure atlas_Mesh__cptr
  module procedure atlas_Mesh__ctor
end interface

!------------------------------------------------------------------------------

