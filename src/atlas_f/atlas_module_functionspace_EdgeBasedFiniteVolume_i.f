! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_functionspace_Nodes) :: atlas_functionspace_EdgeBasedFiniteVolume

! Purpose :
! -------
!   *atlas_functionspace_EdgeBasedFiniteVolume* : Interpretes fields defined in nodes,
!         but in a edge-based finite volume context.
!         The mesh will be adapted to cater for this upon construction

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains

END TYPE atlas_functionspace_EdgeBasedFiniteVolume

interface atlas_functionspace_EdgeBasedFiniteVolume
  module procedure atlas_functionspace_EdgeBasedFiniteVolume__cptr
  module procedure atlas_functionspace_EdgeBasedFiniteVolume__mesh_config
end interface

!------------------------------------------------------------------------------

