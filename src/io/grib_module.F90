! module for testing creation of a Grid
module grib_module
#ifdef HAVE_GRIB
  use grib_api
#endif
  use common_module
  use datastruct_module
  use parallel_module, only: myproc
  implicit none
contains

subroutine write_grib(dstruct,filename)
    implicit none
    type(DataStructure_type), intent(inout), target :: dstruct
    character(len=*), intent(in) :: filename

#ifdef HAVE_GRIB
    integer :: ifile_handle,igrib_handle,iret
    integer :: ihour

    real(kind=jprw), pointer :: loc_depth(:)
    real(kind=jprw), pointer :: loc_velocity(:,:)

    real(kind=jprw), allocatable :: depth(:)
    real(kind=jprw), allocatable :: velocity(:,:)
    real(kind=jprw), allocatable :: Vx(:)
    real(kind=jprw), allocatable :: Vy(:)


    call log_info( "Writing Grib file "//trim(filename) )

    loc_depth => scalar_field_2d("depth",dstruct)
    loc_velocity => vector_field_2d("velocity",dstruct)

    call gather(loc_depth,depth, dstruct)
    call gather(loc_velocity,velocity, dstruct)

    if (myproc .eq. 0) then 

        allocate( Vx( size(depth) ) )
        allocate( Vy( size(depth) ) )

        Vx(:) = velocity(1,:) / depth(:)
        Vy(:) = velocity(2,:) / depth(:)
        ihour = dstruct%time / 3600.

        call grib_open_file(IFILE_HANDLE,'data/gribfield','R')
        call grib_new_from_file(IFILE_HANDLE,IGRIB_HANDLE, IRET)
        call grib_close_file(IFILE_HANDLE)


        call grib_open_file(IFILE_HANDLE,trim(filename),'W')
        call grib_set(IGRIB_HANDLE,'level',1)
        call grib_set(IGRIB_HANDLE,'stepUnits','h')
        call grib_set(IGRIB_HANDLE,'endStep',ihour)

        call grib_set(IGRIB_HANDLE,'paramId',130)
        call grib_set(IGRIB_HANDLE,'values',depth)
        call grib_write(IGRIB_HANDLE,IFILE_HANDLE)

        call grib_set(IGRIB_HANDLE,'paramId',131)
        call grib_set(IGRIB_HANDLE,'values',Vx)
        call grib_write(IGRIB_HANDLE,IFILE_HANDLE)

        call grib_set(IGRIB_HANDLE,'paramId',132)
        call grib_set(IGRIB_HANDLE,'values',Vy)

        call grib_write(IGRIB_HANDLE,IFILE_HANDLE)

        call grib_release(IGRIB_HANDLE)
        call grib_close_file(IFILE_HANDLE)
    end if
#endif
end subroutine write_grib

end module grib_module
