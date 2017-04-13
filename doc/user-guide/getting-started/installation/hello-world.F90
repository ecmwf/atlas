program hello_world

use atlas_module, only : atlas_library, atlas_log

call atlas_library%initialise()
call atlas_log%info("Hello world!")
call atlas_library%finalise()

end program hello_world
