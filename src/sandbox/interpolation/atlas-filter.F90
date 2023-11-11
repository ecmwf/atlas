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

implicit none

contains
    
TYPE(atlas_Redistribution) :: src_dist, tgt_dist
TYPE(atlas_Interpolation)  :: interpolation_st, interpolation_ts
TYPE(atlas_Field)          :: src_field, tgt_field, src_field_tgtpart, tgt_field_srcpart

INTEGER(KIND=JPIM) :: INODE_SIZE
LOGICAL   , POINTER :: LLIS_GHOST(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IPA(:)

PUBLIC :: FILTER_SETUP, FILTER_EXECUTE
PRIVATE :: src_dist, tgt_dist, interpolation_st, interpolation_ts, &
    & src_field, tgt_field, src_field_tgtpart, tgt_field_srcpart

FINAL :: FILTER_FINALISE

CONTAINS

SUBROUTINE FILTER_SETUP(YDGEOMETRY,YDFIELDS,PGMV,PGMVS,PGMVO,PGMVSO)
REAL(JPRB), POINTER :: SRC_VAL(:,:)
TYPE(ATLAS_MESH_NODES) :: NODES
TYPE(ATLAS_FIELD) :: GHOSTFIELD
TYPE(ATLAS_FUNCTIONSPACE_NODECOLUMNS) :: NODES_FS
type(atlas_Grid)           :: src_grid, tgt_grid
type(atlas_Config)         :: interpolation_config
class(atlas_FunctionSpace) :: src_fs, tgt_fs    !! source CellColumns/NodeColumns where source is independent of target

src_grid = atlas_StructuredGrid("O80")
tgt_grid = atlas_StructuredGrid("O40")
src_fs = atlas_functionspace_NodeColumns(grid, atlas_Partitioner("equal_regions"), halo=4)
tgt_fs = atlas_functionspace_NodeColumns(grid, atlas_Partitioner("regular_bands"), halo=2)

! // redistribution setup
src_mesh_tgtpart = altas_Mesh(src_grid, atlas_MatchingPartitioner(tgt_mesh))
tgt_mesh_srcpart = atlas_Mesh(tgt_grid, atlas_MatchingPartitioner(src_mesh))
src_fs_tgtpart = atlas_functionspace_NodeColumns(src_mesh_tgtpart)
tgt_fs_srcpart = atlas_functionspace_NodeColumns(tgt_mesh_srcpart)
src_redist = atlas_Redistribution(src_fs, src_fs_tgtpart)
tgt_redist = atlas_Redistribution(tgt_fs, tgt_fs_srcpart)

! // interpolation setup
interpolation_config = atlas_Config()
call interpolation_config%set("type", "conservative-spherical-polygon")
interpolation_st = atlas_Interpolation(interpolation_config, src_fs_tgtpart, tgt_fs)   !! generate interpolation matrix (divided among the MPI-tasks)
interpolation_ts = atlas_Interpolation(interpolation_config, tgt_fs_srcpart, src_fs)   !! generate interpolation matrix (divided among the MPI-tasks)

! // fields and helper fields
src_field = src_fs%create_field(name="var", kind=atlas_real(JPRB), IGPVARS)
call src_field%data(src_val)

tgt_field = tgt_fs%create_field(name="var", kind=atlas_real(JPRB),IGPVARS)
src_field_tgtpart = src_fs_tgtpart%create_field(name="var_tmp", kind=atlas_real(JPRB),IGPVARS)
tgt_field_srcpart = tgt_fs_srcpart%create_field(name="var_tmp", kind=atlas_real(JPRB),IGPVARS)

! // free memory on all unnecessary Atlas structures by calling 'final' from each object
call src_mesh%final()
call tgt_mesh%final()
call src_mesh_tgtpart%final()
call tgt_mesh_srcpart%final()
call src_fs_tgtpart%final()
call tgt_fs_srcpart%final()
call src_fs%final()
call tgt_fs%final()
CALL GHOSTFIELD%FINAL()
CALL NODES_FS%FINAL()
CALL NODES%FINAL()

!     ------------------------------------------------------------------
END ASSOCIATE

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('FILTER_SETUP',1,ZHOOK_HANDLE)
END SUBROUTINE FILTER_SETUP
end module filter_module

! --------------------------------------------------------------------------------------------

program atlas_filter

use atlas_module
use filter_module

implicit none

type(atlas_Grid) :: grid
type(filter_module) :: filter

grid = atlas_Grid("O40")
print *, "size ", grid%size()

end program atlas_filter
