#ifndef BuildDualMesh_hpp
#define BuildDualMesh_hpp

namespace ecmwf {
  class Mesh;

void build_dual_mesh( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void ecmwf__build_dual_mesh (Mesh* mesh);
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // BuildEdges_hpp
