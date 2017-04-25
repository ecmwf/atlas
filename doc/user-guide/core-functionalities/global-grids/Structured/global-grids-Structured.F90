program main
use atlas_module
implicit none
character(len=1024)         :: string
type(atlas_StructuredGrid)  :: grid

call atlas_library%initialise()

grid = atlas_StructuredGrid( "O32" )

write(string, "(A,I0)") "ny = ", grid%ny()
call atlas_log%info(string)

write(string, "(A,I0)") "nx first = ", grid%nx(1)
call atlas_log%info(string)

write(string, "(A,I0)") "npts = ", grid%size()
call atlas_log%info(string)

write(string, "(A,F8.4)") "y(1)   = ", grid%y(1)
call atlas_log%info(string)

write(string, "(A,F8.4)") "x(1,1) = ", grid%x(1, 1)
call atlas_log%info(string)

call atlas_library%finalise()

end program main
