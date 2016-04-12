program hello_world

use atlas_module, only : &
  & atlas_init, &
  & atlas_finalize, &
  & atlas_log

call atlas_init()
call atlas_log%info("Hello world!")
call atlas_finalize()

end program hello_world
