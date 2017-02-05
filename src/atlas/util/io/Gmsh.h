/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/util/Metadata.h"

namespace atlas { namespace mesh { class Mesh; } }
namespace atlas {
namespace field {
    class Field;
    class FieldSet;
} }

namespace atlas { namespace array { class Array; } }

namespace atlas { namespace util { namespace io { class GmshFortranInterface; } } }
namespace atlas { namespace output { class Gmsh; } }

namespace atlas {
namespace util {
namespace io {

//----------------------------------------------------------------------------------------------------------------------

class Gmsh {

public: // members

  Metadata options;

private: // types

  typedef std::ios_base::openmode openmode;

public: // methods

  Gmsh();

  ~Gmsh();

public:

  mesh::Mesh* read(const eckit::PathName& file_path) const;

  void read(const eckit::PathName& file_path, mesh::Mesh& mesh) const;

  /// Write mesh file
  /// Extra file with available mesh information is written to:
  ///  - filename_info.msh
  void write(const mesh::Mesh& mesh,
             const eckit::PathName& file_path) const;

 /// Write field to file
 ///  Depending on argument "mode", the fields will be appended,
 ///  or existing file will be overwritten
 void write(const field::Field& field,
            const eckit::PathName& file_path,
            openmode mode = std::ios::out) const;

  /// Write fieldset to file using FunctionSpace
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write(const field::FieldSet& fieldset,
             const functionspace::FunctionSpace&,
             const eckit::PathName& file_path,
             openmode mode = std::ios::out) const;

  /// Write field to file using Functionspace
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write(const field::Field& field,
             const functionspace::FunctionSpace&,
             const eckit::PathName& file_path,
             openmode mode = std::ios::out) const;






private:

  /// Write fieldset to file using Nodes functionspace
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write_delegate(const field::FieldSet& fieldset,
                      const functionspace::NodeColumns&,
                      const eckit::PathName& file_path,
                      openmode mode = std::ios::out) const;

  /// Write fieldset to file using StructuredColumns functionspace
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write_delegate(const field::FieldSet& fieldset,
                      const functionspace::StructuredColumns&,
                      const eckit::PathName& file_path,
                      openmode mode = std::ios::out) const;

  /// Write field to file using Nodes functionspace
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write_delegate(const field::Field& field,
             const functionspace::NodeColumns&,
             const eckit::PathName& file_path,
             openmode mode = std::ios::out) const;

  /// Write field to file using StructuredColumns functionspace
  ///  Depending on argument "mode", the fields will be appended,
  ///  or existing file will be overwritten
  void write_delegate(const field::Field& field,
             const functionspace::StructuredColumns&,
             const eckit::PathName& file_path,
             openmode mode = std::ios::out) const;

};

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

extern "C" {
Gmsh* atlas__Gmsh__new();
void atlas__Gmsh__delete(Gmsh* This);
mesh::Mesh* atlas__Gmsh__read(Gmsh* This, char* file_path);
void atlas__Gmsh__write(Gmsh* This, mesh::Mesh* mesh, char* file_path);
mesh::Mesh* atlas__read_gmsh(char* file_path);
void atlas__write_gmsh_mesh(mesh::Mesh* mesh, char* file_path);
void atlas__write_gmsh_fieldset(field::FieldSet* fieldset, functionspace::FunctionSpace* function_space, char* file_path, int mode);
void atlas__write_gmsh_field(field::Field* field, functionspace::FunctionSpace* function_space, char* file_path, int mode);
}

//----------------------------------------------------------------------------------------------------------------------

} // namespace io
} // namespace util
} // namespace atlas

#endif  // Gmsh_h
