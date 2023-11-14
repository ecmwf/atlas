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
!   Filip Vana, Slavko Brdar, Willem Deconinck (Nov 2023)
! 

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------

module filter_module

use atlas_module
use atlas_redistribution_module

implicit none

public :: atlas_Filter, JPRB

private

INTEGER, PARAMETER:: JPRB = SELECTED_REAL_KIND(13,300)

! --------------------------------------------------------------------------------------------
type :: atlas_Filter
private
    type(atlas_FunctionSpace) :: src_fs
    type(atlas_Redistribution):: src_redist, tgt_redist
    type(atlas_Interpolation) :: interpolation_st, interpolation_ts
    type(atlas_Field)         :: tgt_field, src_field_tgtpart, tgt_field_srcpart
contains
    procedure, public :: execute => filter_execute
    procedure, public :: source  => filter_src_fs
    final :: filter_finalise
end type atlas_Filter
! --------------------------------------------------------------------------------------------

interface atlas_Filter
    module procedure atlas_Filter__create
end interface

CONTAINS

! --------------------------------------------------------------------------------------------

function filter_src_fs(this) result(fs)
    class(atlas_Filter), intent(in) :: this
    type(atlas_FunctionSpace) :: fs
    fs = this%src_fs
end function filter_src_fs

! --------------------------------------------------------------------------------------------

function atlas_Filter__create(src_grid, src_mesh) result(this)
    type(atlas_Grid), intent(in) :: src_grid
    type(atlas_Mesh), intent(inout) :: src_mesh
    type(atlas_Filter)           :: this
    type(atlas_Grid)             :: tgt_grid
    type(atlas_MeshGenerator)    :: meshgen
    type(atlas_GridDistribution) :: griddist
    type(atlas_Redistribution)   :: src_redist, tgt_redist
    type(atlas_Mesh)             :: tgt_mesh
    type(atlas_Mesh)             :: src_mesh_tgtpart, tgt_mesh_srcpart
    type(atlas_Config)           :: interpolation_config
    type(atlas_FunctionSpace)    :: src_fs, tgt_fs
    type(atlas_FunctionSpace)    :: src_fs_tgtpart, tgt_fs_srcpart

    tgt_grid = atlas_StructuredGrid("O40")
    meshgen = atlas_MeshGenerator()
    griddist = atlas_GridDistribution(tgt_grid, atlas_Partitioner("regular_bands"))
    tgt_mesh = meshgen%generate(tgt_grid, griddist)
    src_fs = atlas_functionspace_NodeColumns(src_mesh, halo=4)
    tgt_fs = atlas_functionspace_NodeColumns(tgt_mesh, halo=2)
    this%src_fs = src_fs

    ! // redistribution setup
    griddist = atlas_GridDistribution(src_grid, atlas_MatchingPartitioner(tgt_mesh))
    src_mesh_tgtpart = meshgen%generate(src_grid, griddist)
    griddist = atlas_GridDistribution(tgt_grid, atlas_MatchingPartitioner(src_mesh))
    tgt_mesh_srcpart = meshgen%generate(tgt_grid, griddist)
    src_fs_tgtpart = atlas_functionspace_NodeColumns(src_mesh_tgtpart)
    tgt_fs_srcpart = atlas_functionspace_NodeColumns(tgt_mesh_srcpart)

    this%src_redist = atlas_Redistribution(src_fs, src_fs_tgtpart)
    this%tgt_redist = atlas_Redistribution(tgt_fs, tgt_fs_srcpart)

    ! // interpolation setup
    interpolation_config = atlas_Config()
    call interpolation_config%set("type", "conservative-spherical-polygon")
    this%interpolation_st = atlas_Interpolation(interpolation_config, src_fs_tgtpart, tgt_fs)
    this%interpolation_ts = atlas_Interpolation(interpolation_config, tgt_fs_srcpart, src_fs)

    ! // prepare helper fields
    this%tgt_field = tgt_fs%create_field(name="var_tmp", kind=atlas_real(JPRB))
    this%src_field_tgtpart = src_fs_tgtpart%create_field(name="var_tmp", kind=atlas_real(JPRB))
    this%tgt_field_srcpart = tgt_fs_srcpart%create_field(name="var_tmp", kind=atlas_real(JPRB))

    ! // free memory
    call src_mesh%final()
    call tgt_mesh%final()
    call src_mesh_tgtpart%final()
    call tgt_mesh_srcpart%final()
    call src_fs_tgtpart%final()
    call tgt_fs_srcpart%final()
    call src_fs%final()
    call tgt_fs%final()
end function atlas_Filter__create

! --------------------------------------------------------------------------------------------

SUBROUTINE FILTER_EXECUTE(this, src_field)
    class(atlas_Filter), intent(inout) :: this
    TYPE(atlas_Field), intent(inout)  :: src_field

    ASSOCIATE(src_redist=>this%src_redist, src_field_tgtpart=>this%src_field_tgtpart, &
            & tgt_redist=>this%tgt_redist, tgt_field_srcpart=>this%tgt_field_srcpart, &
            & tgt_field=>this%tgt_field)
    ASSOCIATE(interpolation_st=>this%interpolation_st, interpolation_ts=>this%interpolation_ts)

    call src_redist%execute(src_field, src_field_tgtpart)
    call src_field_tgtpart%halo_exchange()
    call interpolation_st%execute(src_field, tgt_field)
    call tgt_redist%execute(tgt_field, tgt_field_srcpart)
    call interpolation_ts%execute(tgt_field_srcpart, src_field)

    END ASSOCIATE
    END ASSOCIATE
END SUBROUTINE FILTER_EXECUTE

! --------------------------------------------------------------------------------------------

SUBROUTINE FILTER_FINALISE(this)
    type(atlas_Filter), intent(inout) :: this

    ASSOCIATE(src_redist=>this%src_redist, src_field_tgtpart=>this%src_field_tgtpart, &
            & tgt_redist=>this%tgt_redist, tgt_field_srcpart=>this%tgt_field_srcpart, &
            & tgt_field=>this%tgt_field)
    ASSOCIATE(interpolation_st=>this%interpolation_st, interpolation_ts=>this%interpolation_ts)

    call src_redist%final()
    call tgt_redist%final()
    call interpolation_st%final()
    call interpolation_ts%final()
    call src_field_tgtpart%final()
    call tgt_field%final()
    call tgt_field_srcpart%final()

    END ASSOCIATE
    END ASSOCIATE
END SUBROUTINE FILTER_FINALISE

! --------------------------------------------------------------------------------------------

end module filter_module

! --------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------

program filtering

use atlas_module
use filter_module, only: JPRB, atlas_Filter

implicit none

    type(atlas_Grid) :: grid
    type(atlas_GridDistribution) :: griddist
    type(atlas_Field) :: sfield
    type(atlas_Filter) :: filter
    type(atlas_FunctionSpace) :: fspace
    type(atlas_Mesh) :: mesh
    type(atlas_MeshGenerator) :: meshgen
    type(atlas_Output) :: gmsh

    real :: start_time, end_time
    real(kind=JPRB), pointer :: sfield_v(:)

    call atlas_library%initialise()

    grid = atlas_StructuredGrid("O80")
    meshgen = atlas_MeshGenerator()
    griddist = atlas_GridDistribution(grid, atlas_Partitioner("equal_regions"))
    mesh = meshgen%generate(grid, griddist)
    gmsh = atlas_output_Gmsh("mesh.msh", "w")
    call gmsh%write(mesh)

    call cpu_time(start_time)
    filter = atlas_Filter(grid, mesh)
    call cpu_time(end_time)
    print *, " filter.setup in seconds: ", start_time - end_time

    fspace = filter%source()
    sfield = fspace%create_field(name="sfield", kind=atlas_real(JPRB))
    call sfield%data(sfield_v)
    sfield_v = 1._JPRB

    call gmsh%write(sfield)

    call cpu_time(start_time)
    call filter.execute(sfield)
    call cpu_time(end_time)
    print *, " filter.exe in seconds: ", start_time - end_time

    call gmsh%write(sfield)

    call atlas_library%finalise()

end program filtering
