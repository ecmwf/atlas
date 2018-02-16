program atlas_acc_fields
use atlas_module
implicit none

type(atlas_Field)  :: field
integer(4), pointer   :: view(:,:,:)
integer :: i,j,k
integer, parameter :: Ni = 2
integer, parameter :: Nj = 3
integer, parameter :: Nk = 5

call atlas_library%initialise()
field = atlas_Field(kind=atlas_integer(4),shape=[Nk,Nj,Ni])
call field%data(view)

!$acc data present(view)
!$acc kernels
do i=1,Ni
  do j=1,Nj
    do k=1,Nk
      view(k,j,i) = i*100*100 + j*100 + k
    enddo
  enddo
enddo
!$acc end kernels
!$acc end data


! We can either use the acc directive, or
! use the field's API to update the host
!     !$acc update host(view)
call field%clone_from_device()


! Note that field%sync_host_device() is not working here...
! Because internal state of host_needs_update has not been changed

do i=1,Ni
  write(0,*) " "
  write(0,*) "i=",i
  write(0,*) view(:,:,i)
enddo

end program
