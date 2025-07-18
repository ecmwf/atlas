! (C) Copyright 2023 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#include "atlas/atlas_f.h"

module atlas_MeshBuilder_module


use fckit_owned_object_module, only: fckit_owned_object

implicit none

private :: fckit_owned_object

public :: atlas_TriangularMeshBuilder

private

!-----------------------------!
! atlas_TriangularMeshBuilder !
!-----------------------------!

!------------------------------------------------------------------------------
TYPE, extends(fckit_owned_object) :: atlas_TriangularMeshBuilder
contains
  procedure, public :: build => atlas_TriangularMeshBuilder__build

#if FCKIT_FINAL_NOT_INHERITING
  final :: atlas_TriangularMeshBuilder__final_auto
#endif

END TYPE atlas_TriangularMeshBuilder

interface atlas_TriangularMeshBuilder
  module procedure atlas_TriangularMeshBuilder__cptr
  module procedure atlas_TriangularMeshBuilder__config
end interface

!------------------------------------------------------------------------------


!========================================================
contains
!========================================================


function atlas_TriangularMeshBuilder__cptr(cptr) result(this)
  use, intrinsic :: iso_c_binding, only: c_ptr
  type(atlas_TriangularMeshBuilder) :: this
  type(c_ptr), intent(in) :: cptr
  call this%reset_c_ptr( cptr )
  call this%return()
end function

function atlas_TriangularMeshBuilder__config(config) result(this)
  use fckit_c_interop_module, only: c_str
  use atlas_MeshBuilder_c_binding
  use atlas_Config_module, only: atlas_Config
  type(atlas_TriangularMeshBuilder) :: this
  type(atlas_Config), intent(in), optional :: config
  call this%reset_c_ptr( atlas__TriangularMeshBuilder__new() )
  call this%return()
end function


!    Mesh operator()(size_t nb_nodes,  const gidx_t node_global_index[], const double x[], const double y[], const double lon[], const double lat[],
!                    size_t nb_triags, const gidx_t triangle_global_index[], const gidx_t triangle_nodes_global_index[],
!                    gidx_t global_index_base) const;

function atlas_TriangularMeshBuilder__build(this, &
    nb_nodes, node_global_index, x, y, lon, lat, &
    nb_triags, triag_global_index, triag_nodes, &
    global_index_base) result(mesh)
  use, intrinsic :: iso_c_binding, only: c_double, c_size_t
  use atlas_MeshBuilder_c_binding
  use atlas_Mesh_module, only: atlas_Mesh
  use atlas_kinds_module, only : ATLAS_KIND_GIDX
  use fckit_array_module, only : array_stride, array_view1d
  type(atlas_Mesh) :: mesh
  class(atlas_TriangularMeshBuilder), intent(in) :: this
  integer, intent(in) :: nb_nodes
  integer(ATLAS_KIND_GIDX), intent(in) :: node_global_index(nb_nodes)
  real(c_double), intent(in), target :: x(:), y(:), lon(:), lat(:)
  integer, intent(in) :: nb_triags
  integer(ATLAS_KIND_GIDX), intent(in), target :: triag_global_index(nb_triags)
  integer(ATLAS_KIND_GIDX), intent(in), target :: triag_nodes(3,nb_triags)
  integer(ATLAS_KIND_GIDX), optional, intent(in) :: global_index_base
  integer(ATLAS_KIND_GIDX) :: global_index_base_
  real(c_double), pointer :: view1d_x(:), view1d_y(:), view1d_lon(:), view1d_lat(:)
  integer(ATLAS_KIND_GIDX), pointer :: view1d_triag_global_index(:), view1d_triag_nodes(:)
  global_index_base_ = 1
  if (present(global_index_base)) then
    global_index_base_ = global_index_base
  endif
  view1d_x => array_view1d(x)
  view1d_y => array_view1d(y)
  view1d_lon => array_view1d(lon)
  view1d_lat => array_view1d(lat)
  view1d_triag_global_index => array_view1d(triag_global_index)
  view1d_triag_nodes => array_view1d(triag_nodes)
  call mesh%reset_c_ptr() ! Somehow needed with PGI/16.7 and build-type "bit"
  mesh = atlas_Mesh( atlas__TriangularMeshBuilder__operator(this%CPTR_PGIBUG_A, &
       & int(nb_nodes,c_size_t), node_global_index, &
       & view1d_x,   view1d_y,   int(array_stride(x,1),c_size_t),   int(array_stride(y,1),c_size_t),   &
       & view1d_lon, view1d_lat, int(array_stride(lon,1),c_size_t), int(array_stride(lat,1),c_size_t), &
       & int(nb_triags,c_size_t), view1d_triag_global_index, view1d_triag_nodes, &
       & global_index_base_ ))
  call mesh%return()
end function


!-------------------------------------------------------------------------------

#if FCKIT_FINAL_NOT_INHERITING
ATLAS_FINAL subroutine atlas_TriangularMeshBuilder__final_auto(this)
  type(atlas_TriangularMeshBuilder), intent(inout) :: this
#if FCKIT_FINAL_DEBUGGING
  write(0,*) "atlas_MeshBuilder__final_auto"
#endif
#if FCKIT_FINAL_NOT_PROPAGATING
  call this%final()
#endif
  FCKIT_SUPPRESS_UNUSED( this )
end subroutine
#endif

! ----------------------------------------------------------------------------------------

end module atlas_MeshBuilder_module
