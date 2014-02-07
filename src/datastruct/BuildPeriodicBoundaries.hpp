#ifndef BuildPeriodicBoundaries_hpp
#define BuildPeriodicBoundaries_hpp
#include <string>
namespace ecmwf {
class Mesh;

void build_periodic_boundaries( Mesh& mesh );

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  void ecmwf__build_periodic_boundaries (Mesh* mesh);
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // BuildPeriodicBoundaries_hpp
