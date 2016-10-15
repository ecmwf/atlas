module atlas_mpi_module
private

public :: atlas_mpi_comm
public :: atlas_mpi_set_comm

contains

function atlas_mpi_comm()
  use fckit_mpi_module
  type(fckit_mpi_comm) :: atlas_mpi_comm
  atlas_mpi_comm = fckit_mpi_comm()
end function atlas_mpi_comm

subroutine atlas_mpi_set_comm(comm)
  use fckit_mpi_module
  integer :: comm
  call fckit_mpi_setCommDefault(comm)
end subroutine atlas_mpi_set_comm

end module atlas_mpi_module
