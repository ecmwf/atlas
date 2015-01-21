
module atlas_mpi_module
#include "atlas/atlas_mpi_fortran.h"
use atlas_mpi_c_binding
public
contains

  subroutine atlas_mpi_init()
    integer :: ierr
    if (.not. atlas_mpi_initialized() ) call MPI_INIT( ierr )
  end subroutine atlas_mpi_init

  subroutine atlas_mpi_finalize()
    integer :: ierr
    logical :: finalized
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
    atlas_mpi_comm = atlas_mpi_Comm_fortran()
  end function atlas_mpi_comm

  subroutine atlas_mpi_set_comm(comm)
    integer :: comm
    call atlas_mpi_Comm_assign(comm)
  end subroutine atlas_mpi_set_comm
  
  subroutine atlas_mpi_barrier(comm)
    integer, optional :: comm
    integer :: ierr
    if( .not. present(comm) ) then
      call MPI_Barrier ( comm, ierr )
    else
      call MPI_Barrier ( atlas_mpi_comm(), ierr )
    endif
  end subroutine atlas_mpi_barrier

  function atlas_mpi_rank(comm)
    integer, optional :: comm
    integer :: ierr
    integer :: MPL_rank
    if( .not. present(comm) ) then
      call MPI_COMM_RANK( comm, MPL_rank, ierr )
    else
      call MPI_COMM_RANK( atlas_mpi_comm(), MPL_rank, ierr )
    endif    
  end function atlas_mpi_rank

  function atlas_mpi_size(comm)
    integer, optional :: comm
    integer :: ierr
    integer :: MPL_size
    if( .not. present(comm) ) then
      call MPI_COMM_SIZE( comm, MPL_size,  ierr )
    else
      call MPI_COMM_SIZE( atlas_mpi_comm(), MPL_size,  ierr )
    endif 
  end function atlas_mpi_size

end module atlas_mpi_module
