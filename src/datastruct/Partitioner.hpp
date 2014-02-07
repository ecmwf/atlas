#ifndef Partitioner_hpp
#define Partitioner_hpp

#include <mpi.h>
#include <vector>
#include <stdexcept>

namespace ecmwf {

class Mesh;

class Partitioner
{
public:
  Partitioner();
  virtual ~Partitioner() {}

public: // methods

  static void partition( Mesh& mesh,
                  int nb_partitions );

private: // methods
private: // data
};

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Partitioner* ecmwf__Partitioner__new ();
  void ecmwf__Partitioner__delete (Partitioner* This);
}
// ------------------------------------------------------------------


} // namespace ecmwf

#endif // Partitioner_hpp
