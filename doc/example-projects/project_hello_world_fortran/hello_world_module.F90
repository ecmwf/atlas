module hello_world_module

use fckit_module
use atlas_module

implicit none

public

contains

subroutine print_hello_world()
  call fckit_log%info("Hello world")
end subroutine

end module
