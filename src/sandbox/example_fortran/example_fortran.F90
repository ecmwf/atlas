
#include "atlas/atlas_f.h"

program example_fortran
use atlas_module
implicit none

integer j
integer :: var(2000)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(j)
do j=1,2000
  var(j) = 1
enddo

call atlas_log%info("Atlas initialized")
call atlas_log%info("version = ["//atlas_version()//"]")

call atlas_library%initialise()

call atlas_log%set_fortran_unit(6)

#if ATLAS_HAVE_OMP
call atlas_log%info("OpenMP enabled")
#endif

call atlas_log%debug("Here is some debugging information")

write(6,'(A)') "exit"
call atlas_library%finalise()
end program
