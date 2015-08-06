! (C) Copyright 2013-2014 ECMWF.


!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_NodesFunctionSpace

! Purpose :
! -------
!   *atlas_NodesFunctionSpace* : Interpretes fields defined in nodes

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------

contains
END TYPE atlas_NodesFunctionSpace

interface atlas_NodesFunctionSpace
  module procedure atlas_NodesFunctionSpace__ctor
end interface

!------------------------------------------------------------------------------
TYPE, extends(object_type) :: atlas_NodesColumnFunctionSpace

! Purpose :
! -------
!   *atlas_NodesColumnFunctionSpace* : Interpretes fields defined in nodes

! Methods :
! -------

! Author :
! ------
!   August-2015 Willem Deconinck     *ECMWF*

!------------------------------------------------------------------------------

contains
END TYPE atlas_NodesColumnFunctionSpace

interface atlas_NodesColumnFunctionSpace
  module procedure atlas_NodesColumnFunctionSpace__ctor
end interface

!------------------------------------------------------------------------------

