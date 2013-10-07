module common_module
implicit none
save
integer, parameter :: jprb = selected_real_kind(13,300)
integer, parameter :: XX = 1
integer, parameter :: YY = 2

contains

function L2norm(array) result(norm)
  real(kind=jprb), intent(in) :: array(:)
  real(kind=jprb) :: norm
  integer :: i
  norm = 0
  do i=1,size(array)
    norm = norm + array(i)**2
  end do
  norm = norm/real(size(array))
end function L2norm

end module common_module