! module for testing creation of a Grid
module grib_module
  use grid_module
  use model_module
  use grib_api
  implicit none
contains

subroutine write_grib_state(state,filename)
    implicit none
    class(State_class), pointer, intent(in) :: state
    character(len=*), intent(in) :: filename
    integer :: ifile_handle,igrib_handle,iret
    integer :: ihour
    class(Field_class), pointer :: field

    real(kind=jprb), dimension(:), allocatable :: field_without_periodic

    integer :: ic, ii, inode, nb_unique_nodes, nb_periodic_nodes, nb_nodes
    integer, dimension(:), allocatable :: unique_nodes

    ! This mess is here because the mesh contains some duplicated points

    nb_nodes = state%function_space%grid%nb_nodes
    nb_periodic_nodes = state%function_space%grid%nb_periodic_nodes
    nb_unique_nodes = nb_nodes - nb_periodic_nodes
    allocate(unique_nodes(nb_unique_nodes))
    ic=0
    ii=1
    do inode=1,nb_nodes
        if(state%function_space%grid%periodic_nodes(ii,2) == inode) then
            if(ii < nb_periodic_nodes) ii=ii+1
        else
            ic=ic+1
            unique_nodes(ic)=inode
        endif
    enddo

    ihour = state%time / 3600.

    call grib_open_file(IFILE_HANDLE,'data/gribfield','R')
    call grib_new_from_file(IFILE_HANDLE,IGRIB_HANDLE, IRET)
    call grib_close_file(IFILE_HANDLE)

    field => state%field("depth")

    allocate( field_without_periodic(nb_unique_nodes) )

    do inode=1, nb_unique_nodes
      field_without_periodic(inode) = field%array( unique_nodes(inode) ,1 )
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
end subroutine write_grib_state

end module grib_module