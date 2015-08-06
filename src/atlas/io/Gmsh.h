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
#include <vector>

#include "atlas/Metadata.h"
#include "atlas/functionspace/NodesFunctionSpace.h"

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

  Mesh* read(const eckit::PathName& file_path) const;

  void read(const eckit::PathName& file_path, Mesh& mesh) const;

  /// Write mesh file
  /// Extra file with available mesh information is written to a different file:
  ///  - filename_info.msh
  void write(const Mesh& mesh, const eckit::PathName& file_path) const;

  /// Write fieldset to file
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write(const FieldSet& fieldset, const functionspace::NodesFunctionSpace&, const eckit::PathName& file_path, openmode mode = std::ios::out) const;

  /// Write field to file
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write(const Field& field, const functionspace::NodesFunctionSpace&, const eckit::PathName& file_path, openmode mode = std::ios::out) const;

public:
  Metadata options;
};

//------------------------------------------------------------------------------------------------------
#define NODESFUNCTIONSPACE functionspace::NodesFunctionSpace
// C wrapper interfaces to C++ routines
extern "C" {
Gmsh* atlas__Gmsh__new();
void atlas__Gmsh__delete(Gmsh* This);
Mesh* atlas__Gmsh__read(Gmsh* This, char* file_path);
void atlas__Gmsh__write(Gmsh* This, Mesh* mesh, char* file_path);
Mesh* atlas__read_gmsh(char* file_path);
void atlas__write_gmsh_mesh(Mesh* mesh, char* file_path);
void atlas__write_gmsh_fieldset(FieldSet* fieldset, NODESFUNCTIONSPACE* function_space, char* file_path, int mode);
void atlas__write_gmsh_field(Field* field, NODESFUNCTIONSPACE* function_space, char* file_path, int mode);
}
#undef NODESFUNCTIONSPACE
//------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas

#endif  // Gmsh_h
