
! ----------------------------------------------------------------------- !
! C-bindings

module atlas_grids_c_binding
implicit none
interface
  subroutine atlas__grids__load() bind(c,name="atlas__grids__load")
  end subroutine
end interface
contains
end module atlas_grids_c_binding

! ----------------------------------------------------------------------- !
! Fortran module

module atlas_grids_module
use atlas_grids_c_binding
contains

subroutine atlas_grids_load()
  call atlas__grids__load()
end subroutine atlas_grids_load

end module atlas_grids_module

