// (C) Copyright 1996-2014 ECMWF.

#ifndef BuildHalo_hpp
#define BuildHalo_hpp
#include <string>
namespace atlas {
class Mesh;

void build_halo( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void atlas__build_halo ( Mesh* mesh, int nb_elems );
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // BuildHalo_hpp
