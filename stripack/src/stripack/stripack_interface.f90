! (C) Copyright 2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module stripack_interface_mod

! iso
use iso_c_binding

use stripack_mod, only: trmesh, trfind, trlist2

implicit none

contains

subroutine stripack_trmesh_c(c_num_nodes, c_xs, c_ys, c_zs, &
    c_list, c_lptr, c_lend, c_lnew, c_near, c_next, c_dist, c_ier) &
           bind(c, name='stripack_trmesh_f90')

integer(c_int), intent(in)    :: c_num_nodes
real(c_double), intent(in)    :: c_xs(c_num_nodes)
real(c_double), intent(in)    :: c_ys(c_num_nodes)
real(c_double), intent(in)    :: c_zs(c_num_nodes)
integer(c_int), intent(out)   :: c_list(6*(c_num_nodes-2))
integer(c_int), intent(out)   :: c_lptr(6*(c_num_nodes-2))
integer(c_int), intent(out)   :: c_lend(c_num_nodes)
integer(c_int), intent(out)   :: c_lnew
integer(c_int), intent(inout) :: c_near(c_num_nodes)
integer(c_int), intent(inout) :: c_next(c_num_nodes)
real(c_double), intent(inout) :: c_dist(c_num_nodes)
integer(c_int), intent(out)   :: c_ier

call trmesh(c_num_nodes, c_xs, c_ys, c_zs, c_list, c_lptr, c_lend, c_lnew, c_near, c_next, c_dist, c_ier)

end subroutine stripack_trmesh_c

subroutine stripack_trfind_c(c_guess_index, c_coords, c_num_nodes, c_xs, c_ys, c_zs, &
    c_list, c_lptr, c_lend, c_bary_coords, c_indices) &
           bind(c, name='stripack_trfind_f90')

integer(c_int), intent(in)    :: c_guess_index
real(c_double), intent(in)    :: c_coords(3)
integer(c_int), intent(in)    :: c_num_nodes
real(c_double), intent(in)    :: c_xs(c_num_nodes)
real(c_double), intent(in)    :: c_ys(c_num_nodes)
real(c_double), intent(in)    :: c_zs(c_num_nodes)
integer(c_int), intent(in)    :: c_list(6*(c_num_nodes-2))
integer(c_int), intent(in)    :: c_lptr(6*(c_num_nodes-2))
integer(c_int), intent(in)    :: c_lend(c_num_nodes)
real(c_double), intent(out)   :: c_bary_coords(3)
integer(c_int), intent(out)   :: c_indices(3)

integer :: guess_index

! convert C 0-based index to Fortran 1-based index
guess_index = c_guess_index + 1

call trfind(guess_index, c_coords, c_num_nodes, c_xs, c_ys, c_zs, &
            c_list, c_lptr, c_lend, &
            c_bary_coords(1), c_bary_coords(2), c_bary_coords(3), &
            c_indices(1), c_indices(2), c_indices(3))

! convert Fortran 1-based indices to C 0-based indices
c_indices = c_indices - 1

end subroutine stripack_trfind_c


subroutine stripack_trlist2_c(c_num_nodes, c_list, c_lptr, c_lend, &
    c_nt, c_ltri, c_ier) &
           bind(c, name='stripack_trlist2_f90')

integer(c_int), intent(in)    :: c_num_nodes
integer(c_int), intent(in)    :: c_list(6*(c_num_nodes-2))
integer(c_int), intent(in)    :: c_lptr(6*(c_num_nodes-2))
integer(c_int), intent(in)    :: c_lend(c_num_nodes)
integer(c_int), intent(out)   :: c_nt
integer(c_int), intent(out)   :: c_ltri(3,2*c_num_nodes-4)
integer(c_int), intent(out)   :: c_ier

integer(c_int) :: jtri, jnode

call trlist2(c_num_nodes,c_list,c_lptr,c_lend,c_nt,c_ltri,c_ier)
do jtri=1,c_nt
    do jnode=1,3
        c_ltri(jnode,jtri) = c_ltri(jnode,jtri) - 1
    enddo
enddo

end subroutine stripack_trlist2_c

end module stripack_interface_mod

