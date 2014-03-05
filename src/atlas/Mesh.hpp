#ifndef Mesh_hpp
#define Mesh_hpp

#include <vector>
#include <map>
#include <string>

namespace atlas {
class FunctionSpace;

class Mesh
{
public:
  virtual ~Mesh();
  // Takes ownership, and will be deleted automatically
  FunctionSpace& add_function_space( FunctionSpace* function_space );
  FunctionSpace& function_space(const std::string& name);
  FunctionSpace& function_space(int idx);
  int nb_function_spaces() { return function_spaces_.size(); }
private:
  std::map< std::string, size_t > index_;
  std::vector< FunctionSpace* > function_spaces_;
};


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Mesh* atlas__Mesh__new ();
  void atlas__Mesh__delete (Mesh* This);
  void atlas__Mesh__add_function_space (Mesh* This, FunctionSpace* function_space); 
  FunctionSpace* atlas__Mesh__function_space (Mesh* This, char* name); 
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // Mesh_hpp
