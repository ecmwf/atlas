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
#include <iostream>

namespace atlas {
class Mesh;
class Field;
class FieldSet;

class Gmsh
{
private:
  typedef std::ios_base::openmode openmode;
public:
  virtual ~Gmsh();

  static Mesh* read(const std::string& file_path);

  static void read(const std::string& file_path, Mesh& mesh );

  /// Write 2 mesh files in ASCII format.
  ///  - filename.msh         will have lon-lat coordinates
  ///  - filename.msh.sphere  will have x-y-z coordinates (sphere)
  /// Extra file with available mesh information is written to a different file:
  ///  - filename_info.msh
  static void write(Mesh& mesh, const std::string& file_path);

  /// Write fieldset to file
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  static void write(FieldSet& fieldset, const std::string& file_path, openmode mode = std::ios::out);

  /// Write field to file
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  static void write(Field& field, const std::string& file_path, openmode mode = std::ios::out);

  /// @todo to be merged with write()

  static void write3dsurf( Mesh& mesh, const std::string& file_path );
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
  void atlas__write_gmsh_mesh (Mesh* mesh, char* file_path);
  void atlas__write_gmsh_fieldset (FieldSet* fieldset, char* file_path, int mode);
  void atlas__write_gmsh_field (Field* field, char* file_path, int mode);
}
// ------------------------------------------------------------------

} // namespace atlas

#endif // Gmsh_hpp
