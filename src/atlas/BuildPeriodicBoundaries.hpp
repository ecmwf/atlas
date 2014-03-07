// (C) Copyright 1996-2014 ECMWF.

#ifndef BuildPeriodicBoundaries_hpp
#define BuildPeriodicBoundaries_hpp
#include <string>
namespace atlas {
class Mesh;

void build_periodic_boundaries( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void atlas__build_periodic_boundaries (Mesh* mesh);
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // BuildPeriodicBoundaries_hpp
