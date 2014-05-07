

#ifndef Partitioner_hpp
#define Partitioner_hpp

#include <vector>
#include <stdexcept>

#include "atlas/atlas_mpi.h"

namespace atlas {

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
  Partitioner* atlas__Partitioner__new ();
  void atlas__Partitioner__delete (Partitioner* This);
}
// ------------------------------------------------------------------


} // namespace atlas

#endif // Partitioner_hpp
