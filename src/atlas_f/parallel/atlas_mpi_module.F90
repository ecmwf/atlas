
module atlas_mpi_module
use atlas_mpi_c_binding
#include "atlas_f/internals/atlas_f_mpi.h"
public
contains

  subroutine atlas_mpi_init()
    integer :: ierr
    if (.not. atlas_mpi_initialized() ) call MPI_INIT( ierr )
  end subroutine atlas_mpi_init

  subroutine atlas_mpi_finalize()
    integer :: ierr
    if (.not. atlas_mpi_finalized()) call MPI_FINALIZE (ierr)
  end subroutine atlas_mpi_finalize

  function atlas_mpi_initialized()
    integer :: ierr
    logical :: atlas_mpi_initialized
    call MPI_INITIALIZED(atlas_mpi_initialized, ierr)
  end function atlas_mpi_initialized

  function atlas_mpi_finalized()
    integer :: ierr
    logical :: atlas_mpi_finalized
    call MPI_FINALIZED(atlas_mpi_finalized, ierr)
  end function atlas_mpi_finalized

  function atlas_mpi_comm()
    integer :: atlas_mpi_comm
    atlas_mpi_comm = atlas_mpi_comm_fortran_communicator()
  end function atlas_mpi_comm

  subroutine atlas_mpi_set_comm(comm)
    integer :: comm
    call atlas_mpi_comm_attach_fortran_communicator(comm)
  end subroutine atlas_mpi_set_comm

  subroutine atlas_mpi_barrier(comm)
    integer, optional :: comm
    integer :: ierr
    if( present(comm) ) then
      call MPI_Barrier ( comm, ierr )
    else
      call MPI_Barrier ( atlas_mpi_comm(), ierr )
    endif
  end subroutine atlas_mpi_barrier

  function atlas_mpi_rank(comm)
    integer, optional :: comm
    integer :: ierr
    integer :: atlas_mpi_rank
    if( present(comm) ) then
      call MPI_COMM_RANK( comm, atlas_mpi_rank, ierr )
    else
      call MPI_COMM_RANK( atlas_mpi_comm(), atlas_mpi_rank, ierr )
    endif
  end function atlas_mpi_rank

  function atlas_mpi_proc(comm)
    integer, optional :: comm
    integer :: ierr
    integer :: irank
    integer :: atlas_mpi_proc
    if( present(comm) ) then
      call MPI_COMM_RANK( comm, irank, ierr )
    else
      call MPI_COMM_RANK( atlas_mpi_comm(), irank, ierr )
    endif
    atlas_mpi_proc = irank + 1
  end function atlas_mpi_proc

  function atlas_mpi_size(comm)
    integer, optional :: comm
    integer :: ierr
    integer :: atlas_mpi_size
    if( present(comm) ) then
      call MPI_COMM_SIZE( comm, atlas_mpi_size,  ierr )
    else
      call MPI_COMM_SIZE( atlas_mpi_comm(), atlas_mpi_size,  ierr )
    endif
  end function atlas_mpi_size

end module atlas_mpi_module
