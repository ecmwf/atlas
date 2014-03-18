// (C) Copyright 1996-2014 ECMWF.

#ifndef BuildDualMesh_hpp
#define BuildDualMesh_hpp

namespace atlas {
  class Mesh;

void build_dual_mesh( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void atlas__build_dual_mesh (Mesh* mesh);
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // BuildEdges_hpp
