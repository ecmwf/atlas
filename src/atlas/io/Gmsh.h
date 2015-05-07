/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef Gmsh_h
#define Gmsh_h

#include <string>
#include <iostream>

#include "atlas/Metadata.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

class Mesh;
class Field;
class FieldSet;

namespace io {

//------------------------------------------------------------------------------------------------------

class Gmsh {
private:

  typedef std::ios_base::openmode openmode;

public:

  Gmsh();

  virtual ~Gmsh();

  static Mesh* read(const std::string& file_path);

  static void read(const std::string& file_path, Mesh& mesh );

  /// Write mesh file
  /// Extra file with available mesh information is written to a different file:
  ///  - filename_info.msh
  void write(Mesh& mesh, const std::string& file_path) const;

  /// Write fieldset to file
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write(FieldSet& fieldset, const std::string& file_path, openmode mode = std::ios::out) const;

  /// Write field to file
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write(Field& field, const std::string& file_path, openmode mode = std::ios::out) const;

  /// @todo to be merged with write()
  static void write3dsurf(const Mesh& mesh, const std::string& file_path );

public:

  Metadata options;

public: // this should really belong in options

  std::vector<long> levels;

};

//------------------------------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------------------------------

} // namespace io
} // namespace atlas

#endif // Gmsh_h
