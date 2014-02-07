#ifndef BuildHalo_hpp
#define BuildHalo_hpp
#include <string>
namespace ecmwf {
class Mesh;

void build_halo( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void ecmwf__build_halo ( Mesh* mesh, int nb_elems );
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // BuildHalo_hpp
