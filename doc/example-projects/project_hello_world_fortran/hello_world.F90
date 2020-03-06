program main

use hello_world_module

implicit none

call atlas_library%initialise()

call print_hello_world()

call atlas_library%finalise()

end program
