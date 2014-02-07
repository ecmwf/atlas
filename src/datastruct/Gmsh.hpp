#ifndef Gmsh_hpp
#define Gmsh_hpp
#include <string>
namespace ecmwf {
class Mesh;

class Gmsh
{
public:
  virtual ~Gmsh();

  static Mesh& read(const std::string& file_path);

  static void write(Mesh& mesh, const std::string& file_path);

private:
  enum { LINE=1, TRIAG=2, QUAD=3, POINT=15 };
};


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Gmsh* ecmwf__Gmsh__new ();
  void ecmwf__Gmsh__delete (Gmsh* This);
  Mesh* ecmwf__Gmsh__read (Gmsh* This, char* file_path);
  void ecmwf__Gmsh__write (Gmsh* This, Mesh* mesh, char* file_path);
  Mesh* ecmwf__read_gmsh (char* file_path);
}
// ------------------------------------------------------------------

} // namespace ecmwf

#endif // Gmsh_hpp
