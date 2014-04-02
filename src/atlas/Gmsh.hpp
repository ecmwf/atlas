// (C) Copyright 1996-2014 ECMWF.

#ifndef Gmsh_hpp
#define Gmsh_hpp

#include <string>

namespace atlas {
class Mesh;

class Gmsh
{
public:
  virtual ~Gmsh();

    /// @warning should return pointer, not reference so it can be deleted
    static Mesh* read(const std::string& file_path);

    static void write(Mesh& mesh, const std::string& file_path);

    /// @todo to be merged with write()

    static void write3dsurf( Mesh& mesh, const std::string& file_path );

private:
  enum { LINE=1, TRIAG=2, QUAD=3, POINT=15 };
};


// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines
extern "C" 
{
  Gmsh* atlas__Gmsh__new ();
  void atlas__Gmsh__delete (Gmsh* This);
  Mesh* atlas__Gmsh__read (Gmsh* This, char* file_path);
  void atlas__Gmsh__write (Gmsh* This, Mesh* mesh, char* file_path);
  Mesh* atlas__read_gmsh (char* file_path);
  void atlas__write_gmsh (Mesh* mesh, char* file_path);
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // Gmsh_hpp
