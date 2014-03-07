// (C) Copyright 1996-2014 ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation nor
// does it submit to any jurisdiction.


#ifndef Partitioner_hpp
#define Partitioner_hpp

#include <mpi.h>
#include <vector>
#include <stdexcept>

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
