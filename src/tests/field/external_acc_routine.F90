! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

subroutine external_acc_routine(view)

  implicit none
  real(4), pointer, intent(inout) :: view(:,:)

  !$acc kernels present(view)
  view(1,1) = 5.
  !$acc end kernels

end subroutine external_acc_routine


subroutine external_acc_routine_devptr(dview)

  implicit none
  real(4), pointer, intent(inout) :: dview(:,:)

  !$acc kernels deviceptr(dview)
  dview(1,1) = 6.
  !$acc end kernels

end subroutine external_acc_routine_devptr
