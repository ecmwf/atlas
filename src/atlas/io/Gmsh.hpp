/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */



#ifndef Gmsh_hpp
#define Gmsh_hpp

#include <string>

namespace atlas {
class Mesh;

class Gmsh
{
public:
  virtual ~Gmsh();

    static Mesh* read(const std::string& file_path);

    static void read(const std::string& file_path, Mesh& mesh );

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
