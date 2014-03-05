#ifndef BuildEdges_hpp
#define BuildEdges_hpp
#include <string>
namespace atlas {
class Mesh;

void build_edges( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void atlas__build_edges (Mesh* mesh);
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // BuildEdges_hpp
