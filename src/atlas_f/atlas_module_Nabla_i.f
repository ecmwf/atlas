! (C) Copyright 2013-2015 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(atlas_RefCounted) :: atlas_Nabla

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
  procedure, public :: delete => atlas_Nabla__delete
  procedure, public :: gradient => atlas_nabla__gradient
  procedure, public :: divergence => atlas_nabla__divergence
  procedure, public :: curl => atlas_nabla__curl
  procedure, public :: laplacian => atlas_nabla__laplacian

END TYPE atlas_Nabla

interface atlas_Nabla
  module procedure atlas_Nabla__cptr
  module procedure atlas_Nabla__functionspace_config
end interface

!------------------------------------------------------------------------------
