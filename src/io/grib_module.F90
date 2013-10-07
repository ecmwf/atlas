! module for testing creation of a Grid
module grib_module
#ifdef ENABLE_GRIB
  use grib_api
#endif
  use common_module
  use datastruct_module
  implicit none
contains

subroutine write_grib(geom,filename)
    implicit none
    type(Geometry_class), intent(inout), target :: geom
    character(len=*), intent(in) :: filename

#ifdef ENABLE_GRIB
    integer :: ifile_handle,igrib_handle,iret
    integer :: ihour

    real(kind=jprb), dimension(:), allocatable :: field_without_periodic
    real(kind=jprb), dimension(:), pointer :: field

    integer :: ic, ii, inode, nb_unique_nodes, nb_ghost_nodes, nb_nodes
    integer, dimension(:), allocatable :: unique_nodes

    ! This mess is here because the mesh contains some duplicated points

    nb_nodes = geom%nb_nodes
    nb_ghost_nodes = geom%nb_ghost_nodes
    nb_unique_nodes = nb_nodes - nb_ghost_nodes
    allocate(unique_nodes(nb_unique_nodes))
    ic=0
    ii=1
    do inode=1,nb_nodes
        if(geom%ghost_nodes(ii,2) == inode) then
            if(ii < nb_ghost_nodes) ii=ii+1
        else
            ic=ic+1
            unique_nodes(ic)=inode
        endif
    enddo

    ihour = geom%fields%time / 3600.

    call grib_open_file(IFILE_HANDLE,'data/gribfield','R')
    call grib_new_from_file(IFILE_HANDLE,IGRIB_HANDLE, IRET)
    call grib_close_file(IFILE_HANDLE)

    field => scalar_field("depth",geom)

    allocate( field_without_periodic(nb_unique_nodes) )

    do inode=1, nb_unique_nodes
      field_without_periodic(inode) = field( unique_nodes(inode) )
    end do

    call grib_open_file(IFILE_HANDLE,trim(filename),'W')
    call grib_set(IGRIB_HANDLE,'level',1)
    call grib_set(IGRIB_HANDLE,'paramId',130)
    call grib_set(IGRIB_HANDLE,'stepUnits','h')
    call grib_set(IGRIB_HANDLE,'endStep',ihour)
    call grib_set(IGRIB_HANDLE,'values',field_without_periodic)
    call grib_write(IGRIB_HANDLE,IFILE_HANDLE)

    call grib_release(IGRIB_HANDLE)
    call grib_close_file(IFILE_HANDLE)

#else
    write(0,*) "Grib format is not enabled, compile with '-DENABLE_GRIB'"
#endif
end subroutine write_grib

end module grib_module