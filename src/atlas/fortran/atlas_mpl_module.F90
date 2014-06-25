
module atlas_mpl_module
#include "atlas/atlas_mpi_fortran.h"

contains

  subroutine MPL_init()
    integer :: ierr
    logical :: initialized
    call MPI_INITIALIZED(initialized, ierr)
    if (.not. initialized) call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myproc, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc,  ierr )
  end subroutine MPL_init

  subroutine MPL_finalize()
    integer :: ierr
    logical :: finalized
    call MPI_FINALIZED(finalized, ierr)
    if (.not. finalized) call MPI_FINALIZE (ierr)
  end subroutine MPL_finalize

  subroutine MPL_barrier()
    integer :: ierr
    call MPI_Barrier ( MPI_COMM_WORLD, ierr );
  end subroutine

end module atlas_mpl_module