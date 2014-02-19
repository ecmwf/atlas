#ifndef MPL_hpp
#define MPL_hpp

#include <mpi.h>

namespace ecmwf {

namespace MPL {

  template<typename DATA_TYPE>
  MPI_Datatype TYPE();
  template<> inline MPI_Datatype TYPE<int>()    { return MPI_INT; }
  template<> inline MPI_Datatype TYPE<float>()  { return MPI_FLOAT; }
  template<> inline MPI_Datatype TYPE<double>() { return MPI_DOUBLE; }

  inline void init(int argc=0, char *argv[]=0)
  {
    int ierr;
    int initialized;
    ierr = MPI_Initialized( &initialized );
    if( !initialized )
      ierr = MPI_Init(&argc,&argv);
  }

  inline void finalize()
  {
    int ierr;
    int finalized;
    ierr = MPI_Finalized( &finalized );
    if( !finalized )
      ierr = MPI_Finalize();
  }

  inline int rank()
  {
    int rank;
    int ierr = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    return rank;
  }

  inline int size()
  {
    int nproc;
    int ierr = MPI_Comm_size( MPI_COMM_WORLD, &nproc );
    return nproc;
  }

} // namespace MPL

} // namepsace ecmwf

#endif // MPL_hpp
