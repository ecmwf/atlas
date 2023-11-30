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
!   Filip Vana, Slavko Brdar, Willem Deconinck (ECMWF, Nov 2023)
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
    type(atlas_Mesh)           :: src_mesh
    type(atlas_FunctionSpace)  :: src_fs
    type(atlas_Redistribution) :: src_redist, tgt_redist
    type(atlas_Interpolation)  :: interpolation_st, interpolation_ts
    type(atlas_Field)          :: tgt_field, src_field_tgtpart, tgt_field_srcpart
    character(len=64)          :: data_loc
contains
    procedure, public :: execute => filter_execute
    procedure, public :: fspace  => filter_src_fs
    procedure, public :: mesh    => filter_src_mesh
    final :: filter_finalise
end type atlas_Filter
! --------------------------------------------------------------------------------------------

interface atlas_Filter
    module procedure atlas_Filter__create
end interface

CONTAINS

! --------------------------------------------------------------------------------------------

function filter_src_fs(this) result(fs)
    class(atlas_Filter), intent(inout) :: this
    type(atlas_FunctionSpace) :: fs
    fs = this%src_fs
end function filter_src_fs

! --------------------------------------------------------------------------------------------

function filter_src_mesh(this) result(mesh)
    class(atlas_Filter), intent(inout) :: this
    type(atlas_Mesh) :: mesh
    mesh = this%src_mesh
end function filter_src_mesh

! --------------------------------------------------------------------------------------------

function atlas_Filter__create(src_grid, partitioner, tgt_grid_name, data_loc_in) result(this)
    class(atlas_Grid), intent(in)          :: src_grid
    character(len=*), intent(in)          :: tgt_grid_name
    type(atlas_Partitioner), intent(in)    :: partitioner
    character(len=*), intent(in), optional :: data_loc_in
    type(atlas_Filter)           :: this

    type(atlas_Grid)             :: tgt_grid
    type(atlas_Partitioner)      :: mpart
    type(atlas_MeshGenerator)    :: meshgen
    type(atlas_GridDistribution) :: griddist
    type(atlas_Redistribution)   :: src_redist, tgt_redist
    type(atlas_Mesh)             :: src_mesh, tgt_mesh
    type(atlas_Mesh)             :: src_mesh_tgtpart, tgt_mesh_srcpart
    type(atlas_Config)           :: interpolation_config
    type(atlas_FunctionSpace)    :: src_fs, tgt_fs
    type(atlas_FunctionSpace)    :: src_fs_tgtpart, tgt_fs_srcpart

    this%data_loc = "CellColumns"
    if (present(data_loc_in)) this%data_loc = trim(data_loc_in)
    tgt_grid = atlas_StructuredGrid(tgt_grid_name)
    meshgen = atlas_MeshGenerator()
    griddist = atlas_GridDistribution(src_grid, partitioner)
    src_mesh = meshgen%generate(src_grid, griddist)
    griddist = atlas_GridDistribution(tgt_grid, partitioner)
    tgt_mesh = meshgen%generate(tgt_grid, griddist)
    if (trim(this%data_loc) == "CellColumns") then
        src_fs = atlas_functionspace_CellColumns(src_mesh, halo=0)
        tgt_fs = atlas_functionspace_CellColumns(tgt_mesh, halo=0)
    else if (trim(this%data_loc) == "NodeColumns") then
        ! NodeColumns require halo>=1 for ConservativeSphericalInterpolation
        src_fs = atlas_functionspace_NodeColumns(src_mesh, halo=1)
        tgt_fs = atlas_functionspace_NodeColumns(tgt_mesh, halo=1)
    else if (trim(this%data_loc) == "StructuredColumns") then
        ! NOTE: ConservativeSphericalInterpolation does not support StructuredColumns
        !       please use: call interpolation_config%set("type", "nearest-neighbour")
        !       or switch to functionspace_NodeColumns instead
        src_fs = atlas_functionspace_StructuredColumns(src_grid, partitioner, halo=1)
        tgt_fs = atlas_functionspace_StructuredColumns(tgt_grid, partitioner, halo=1)
    else
        stop "Unknown function space: "//data_loc_in
    end if
    this%src_mesh = src_mesh
    this%src_fs = src_fs

    ! // redistribution setup
    mpart = atlas_MatchingPartitioner(tgt_mesh)
    griddist = atlas_GridDistribution(src_grid, mpart)
    src_mesh_tgtpart = meshgen%generate(src_grid, griddist)
    mpart = atlas_MatchingPartitioner(src_mesh)
    griddist = atlas_GridDistribution(tgt_grid, mpart)
    tgt_mesh_srcpart = meshgen%generate(tgt_grid, griddist)
    if (trim(this%data_loc) == "CellColumns") then
        src_fs_tgtpart = atlas_functionspace_CellColumns(src_mesh_tgtpart, halo=2)
        tgt_fs_srcpart = atlas_functionspace_CellColumns(tgt_mesh_srcpart, halo=1)
    else if (trim(this%data_loc) == "NodeColumns") then
        ! NOTE: sources have to cover the added half-cell size layer on target partitions
        !         For the backward remapping, target cell have to covert the added half-cell
        !         size layer on the source partitions
        !       In praxis, O80 needs 2 cell-layers to cover 1 cell-layer of O40, so O80-halo should be 2 + 1
        !         O40 needs 1 cell-layers to cover 1 cell-layer of O80. O40-halo should be 1 + 1
        src_fs_tgtpart = atlas_functionspace_NodeColumns(src_mesh_tgtpart, halo=9)
        tgt_fs_srcpart = atlas_functionspace_NodeColumns(tgt_mesh_srcpart, halo=2)
    else if (trim(this%data_loc) == "StructuredColumns") then
        src_fs_tgtpart = atlas_functionspace_StructuredColumns(src_grid, partitioner, halo=2)
        tgt_fs_srcpart = atlas_functionspace_StructuredColumns(tgt_grid, partitioner, halo=1)
    endif

    this%src_redist = atlas_Redistribution(src_fs, src_fs_tgtpart)
    this%tgt_redist = atlas_Redistribution(tgt_fs, tgt_fs_srcpart)

    ! // interpolation setup
    interpolation_config = atlas_Config()
    call interpolation_config%set("type", "conservative-spherical-polygon")
    this%interpolation_st = atlas_Interpolation(interpolation_config, src_fs_tgtpart, tgt_fs)
    this%interpolation_ts = atlas_Interpolation(interpolation_config, tgt_fs_srcpart, src_fs)

    ! // prepare helper fields
    this%tgt_field = tgt_fs%create_field(kind=atlas_real(JPRB))
    this%src_field_tgtpart = src_fs_tgtpart%create_field(kind=atlas_real(JPRB))
    this%tgt_field_srcpart = tgt_fs_srcpart%create_field(kind=atlas_real(JPRB))

    ! // free memory
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
    call interpolation_st%execute(src_field_tgtpart, tgt_field)
    call tgt_redist%execute(tgt_field, tgt_field_srcpart)
    call interpolation_ts%execute(tgt_field_srcpart, src_field)

    END ASSOCIATE
    END ASSOCIATE
END SUBROUTINE FILTER_EXECUTE

! --------------------------------------------------------------------------------------------

SUBROUTINE FILTER_FINALISE(this)
    type(atlas_Filter), intent(inout) :: this

    ASSOCIATE(src_mesh=>this%src_mesh, src_fs => this%src_fs)
    ASSOCIATE(src_redist=>this%src_redist, src_field_tgtpart=>this%src_field_tgtpart, &
            & tgt_redist=>this%tgt_redist, tgt_field_srcpart=>this%tgt_field_srcpart, &
            & tgt_field=>this%tgt_field)
    ASSOCIATE(interpolation_st=>this%interpolation_st, interpolation_ts=>this%interpolation_ts)

    ! destroy atlas structures created by this module
    call src_mesh%final()
    call src_fs%final()
    call src_redist%final()
    call tgt_redist%final()
    call interpolation_st%final()
    call interpolation_ts%final()
    call src_field_tgtpart%final()
    call tgt_field%final()
    call tgt_field_srcpart%final()

    END ASSOCIATE
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

use fckit_mpi_module
use, intrinsic :: iso_c_binding

implicit none

    type(atlas_StructuredGrid)   :: grid
    type(atlas_Partitioner)      :: partitioner
    type(atlas_GridDistribution) :: griddist
    type(atlas_Field)            :: sfield, field_lonlat
    type(atlas_Filter)           :: filter
    type(atlas_FunctionSpace)    :: fspace
    type(atlas_Mesh)             :: mesh
    type(atlas_MeshGenerator)    :: meshgen
    type(atlas_Output)           :: gmsh
    type(atlas_mesh_Nodes)       :: nodes

    real                         :: start_time, end_time, mintime, maxtime
    real(kind=JPRB), pointer     :: sfield_v(:), lonlat(:,:)
    integer                      :: inode, nb_nodes

    type(fckit_mpi_comm) :: comm

    call atlas_library%initialise()
    comm = fckit_mpi_comm()

    grid = atlas_StructuredGrid("O320")
    partitioner = atlas_Partitioner("equal_regions")

    call cpu_time(start_time)
    filter = atlas_Filter(grid, partitioner, "O40", "NodeColumns") ! TODO: CellColumns has a problem
    call cpu_time(end_time)
    maxtime = end_time - start_time
    mintime = maxtime
    call comm%allreduce(mintime, fckit_mpi_min())
    call comm%allreduce(maxtime, fckit_mpi_max())
    if(comm%rank()==0) print *, "  filter.setup in seconds (min, max): ", mintime, maxtime

    mesh = filter%mesh()
    gmsh = atlas_output_Gmsh("mesh.msh", "w")
    call gmsh%write(mesh)

    fspace = filter%fspace()
    sfield = fspace%create_field(name="unfiltered", kind=atlas_real(JPRB))

    ! initial data
    call sfield%data(sfield_v)
    nodes = mesh%nodes()
    nb_nodes = nodes%size()
    field_lonlat = nodes%lonlat()
    call field_lonlat%data(lonlat)
    sfield_v(:) = MDPI_gulfstream(lonlat(1,:), lonlat(2,:))
    call sfield%halo_exchange()
    call gmsh%write(sfield)

    call cpu_time(start_time)
    call filter%execute(sfield)
    call cpu_time(end_time)
    maxtime = end_time - start_time
    mintime = maxtime
    call comm%allreduce(mintime, fckit_mpi_min())
    call comm%allreduce(maxtime, fckit_mpi_max())
    if(comm%rank()==0) print *, "  filter.exe in seconds (min, max): ", mintime, maxtime

    call sfield%rename("filtered")
    call sfield%halo_exchange()
    call gmsh%write(sfield)

    ! output the difference of the unfiltered and the filtered field
    sfield_v(:) = sfield_v(:) - MDPI_gulfstream(lonlat(1,:), lonlat(2,:))
    call sfield%rename("filt-unfilt")
    call sfield%halo_exchange()
    call gmsh%write(sfield)
    if(comm%rank()==0) print *, "  output written to: mesh.msh (use GMesh: gmsh mesh.msh)"

    call gmsh%final()
    call nodes%final()
    call partitioner%final()
    call grid%final()
    call atlas_library%finalise()
end program filtering
