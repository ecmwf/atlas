#ifndef mesh_hpp
#define mesh_hpp

#include <vector>
#include <map>
#include <string>

namespace ecmwf {
class FunctionSpace;

class Mesh
{
public:
  virtual ~Mesh();
  // Takes ownership, and will be deleted automatically
  void add_function_space( FunctionSpace* function_space );
  FunctionSpace& function_space(const std::string& name);
private:
  std::map< std::string, size_t > index_;
  std::vector< FunctionSpace* > function_spaces_;
};


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Mesh* ecmwf__Mesh__new ();
  void ecmwf__Mesh__delete (Mesh* This);
  void ecmwf__Mesh__add_function_space (Mesh* This, FunctionSpace* function_space); 
  FunctionSpace* ecmwf__Mesh__function_space (Mesh* This, char* name); 
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // Mesh_hpp
