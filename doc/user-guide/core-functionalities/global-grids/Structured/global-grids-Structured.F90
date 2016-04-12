program main

use atlas_module

character(len=1024)     :: string
character(len=1024)     :: gridID
type(atlas_grid_Structured) :: grid
call atlas_init()

call atlas_resource("--grid", "N32", gridID)

grid = atlas_grid_Structured(gridID)

write(string, "(A,I0)") "nlat = ", grid%nlat()
call atlas_log%info(string)

write(string, "(A,I0)") "nlon = ", grid%nlon(1)
call atlas_log%info(string)

write(string, "(A,I0)") "npts = ", grid%npts()
call atlas_log%info(string)

write(string, "(A,F8.4)") "lat(1)   = ", grid%lat(1)
call atlas_log%info(string)

write(string, "(A,F8.4)") "lon(1,1) = ", grid%lon(1, 1)
call atlas_log%info(string)

call atlas_finalize()

end program main
