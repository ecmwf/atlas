#ifndef BuildEdges_hpp
#define BuildEdges_hpp
#include <string>
namespace ecmwf {
class Mesh;

void build_edges( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void ecmwf__build_edges ( Mesh* mesh );
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // BuildEdges_hpp
