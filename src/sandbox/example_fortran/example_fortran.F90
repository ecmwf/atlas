
program example_fortran 
use atlas_module
implicit none
 
call atlas_init()
call logger%info("Check OK")
call atlas_finalize()

end program
