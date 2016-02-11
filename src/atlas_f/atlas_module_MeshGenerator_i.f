! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_MeshGenerator

! Purpose :
! -------
!   *Nabla* :

! Methods :
! -------

! Author :
! ------
!   October-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------
contains
  procedure, public :: delete => atlas_MeshGenerator__delete
  procedure, public :: copy => atlas_MeshGenerator__copy
  procedure, public :: generate => atlas_MeshGenerator__generate

END TYPE atlas_MeshGenerator

interface atlas_MeshGenerator
  module procedure atlas_MeshGenerator__cptr
  module procedure atlas_MeshGenerator__name_config
end interface

interface atlas_ReducedGridMeshGenerator
  module procedure atlas_ReducedGridMeshGenerator__config
end interface

!------------------------------------------------------------------------------
