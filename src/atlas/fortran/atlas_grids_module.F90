
module atlas_grids_module

use atlas_grids_c_binding

contains

subroutine atlas_grids_load()
  call atlas__grids__load()
end subroutine atlas_grids_load

end module atlas_grids_module

