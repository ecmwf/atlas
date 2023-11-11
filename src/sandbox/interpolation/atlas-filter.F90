!
! (C) Copyright 1996- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
 

!
! Purpose:
!   Smooth GMV and GMVS fields by conservative remapping to a lower resolution grid (and back)
! 
! Authors:
!   Filip Vana, Slavko Brdar, Willem Deconinck
! 

! --------------------------------------------------------------------------------------------

module filter_module

use atlas_module
use atlas_redistribution_module

implicit none

INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300)

TYPE(atlas_Redistribution) :: src_redist, tgt_redist
TYPE(atlas_Interpolation)  :: interpolation_st, interpolation_ts
TYPE(atlas_Field)          :: tgt_field, src_field_tgtpart, tgt_field_srcpart

PUBLIC :: FILTER_SETUP, FILTER_EXECUTE, FILTER_FINALISE
PRIVATE :: src_redist, tgt_redist, interpolation_st, interpolation_ts, &
    & src_field_tgtpart, tgt_field_srcpart

CONTAINS

! --------------------------------------------------------------------------------------------

SUBROUTINE FILTER_SETUP()
    type(atlas_Grid)           :: src_grid, tgt_grid
    type(atlas_MeshGenerator)  :: meshgen
    type(atlas_GridDistribution) :: griddist
    type(atlas_Mesh)           :: src_mesh, tgt_mesh
    type(atlas_Mesh)           :: src_mesh_tgtpart, tgt_mesh_srcpart
    type(atlas_Config)         :: interpolation_config
    type(atlas_FunctionSpace) :: src_fs, tgt_fs
    type(atlas_FunctionSpace) :: src_fs_tgtpart, tgt_fs_srcpart
    type(atlas_Redistribution) :: src_redist, tgt_redist
    real(kind=JPRB), pointer :: src_val(:)

    src_grid = atlas_StructuredGrid("O80")
    tgt_grid = atlas_StructuredGrid("O40")
    meshgen = atlas_MeshGenerator()
    griddist = atlas_GridDistribution(src_grid, atlas_Partitioner("equal_regions"))
    src_mesh = meshgen%generate(src_grid, griddist)
    griddist = atlas_GridDistribution(tgt_grid, atlas_Partitioner("regular_bands"))
    tgt_mesh = meshgen%generate(tgt_grid, griddist)
    src_fs = atlas_functionspace_NodeColumns(src_mesh, halo=4)
    tgt_fs = atlas_functionspace_NodeColumns(tgt_mesh, halo=2)

    ! // redistribution setup
    griddist = atlas_GridDistribution(src_grid, atlas_MatchingPartitioner(tgt_mesh))
    src_mesh_tgtpart = meshgen%generate(src_grid, griddist)
    griddist = atlas_GridDistribution(tgt_grid, atlas_MatchingPartitioner(src_mesh))
    tgt_mesh_srcpart = meshgen%generate(tgt_grid, griddist)
    src_fs_tgtpart = atlas_functionspace_NodeColumns(src_mesh_tgtpart)
    tgt_fs_srcpart = atlas_functionspace_NodeColumns(tgt_mesh_srcpart)
    src_redist = atlas_Redistribution(src_fs, src_fs_tgtpart)
    tgt_redist = atlas_Redistribution(tgt_fs, tgt_fs_srcpart)

    ! // interpolation setup
    interpolation_config = atlas_Config()
    call interpolation_config%set("type", "conservative-spherical-polygon")
    interpolation_st = atlas_Interpolation(interpolation_config, src_fs_tgtpart, tgt_fs)
    interpolation_ts = atlas_Interpolation(interpolation_config, tgt_fs_srcpart, src_fs)

    ! // prepare helper fields
    tgt_field = tgt_fs%create_field(name="var_tmp", kind=atlas_real(JPRB))
    src_field_tgtpart = src_fs_tgtpart%create_field(name="var_tmp", kind=atlas_real(JPRB))
    tgt_field_srcpart = tgt_fs_srcpart%create_field(name="var_tmp", kind=atlas_real(JPRB))

    ! // free memory
    call src_mesh%final()
    call tgt_mesh%final()
    call src_mesh_tgtpart%final()
    call tgt_mesh_srcpart%final()
    call src_fs_tgtpart%final()
    call tgt_fs_srcpart%final()
    call src_fs%final()
    call tgt_fs%final()
END SUBROUTINE FILTER_SETUP

! --------------------------------------------------------------------------------------------

SUBROUTINE FILTER_EXECUTE(src_field)
    TYPE(atlas_Field), intent(inout)  :: src_field

    call src_redist%execute(src_field, src_field_tgtpart)
    call src_field_tgtpart%halo_exchange()
    call interpolation_st%execute(src_field, tgt_field)
    call tgt_redist%execute(tgt_field, tgt_field_srcpart)
    call interpolation_ts%execute(tgt_field_srcpart, src_field)
END SUBROUTINE FILTER_EXECUTE

! --------------------------------------------------------------------------------------------

SUBROUTINE FILTER_FINALISE()
    call src_redist%final()
    call tgt_redist%final()
    call interpolation_st%final()
    call interpolation_ts%final()
    call src_field_tgtpart%final()
    call tgt_field%final()
    call tgt_field_srcpart%final()
END SUBROUTINE FILTER_FINALISE

! --------------------------------------------------------------------------------------------

end module filter_module

! --------------------------------------------------------------------------------------------

program atlas_filter

use atlas_module
use filter_module

implicit none

type(atlas_Grid) :: grid

call filter_setup()
grid = atlas_Grid("O40")
print *, "size ", grid%size()

end program atlas_filter
