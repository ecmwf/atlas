/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/util/Metadata.h"

namespace atlas {
class Field;
class FieldSet;
class Mesh;
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field
namespace functionspace {
class FunctionSpaceImpl;
}
}  // namespace atlas

namespace atlas {
namespace array {
class Array;
}
}  // namespace atlas

namespace atlas {
namespace util {
namespace io {
class GmshFortranInterface;
}
}  // namespace util
}  // namespace atlas
namespace atlas {
namespace output {
class Gmsh;
}
}  // namespace atlas

namespace atlas {
namespace output {
namespace detail {

//----------------------------------------------------------------------------------------------------------------------

class GmshIO {
public:  // members
    using Metadata = util::Metadata;
    Metadata options;

private:  // types
    typedef std::ios_base::openmode openmode;

public:  // methods
    GmshIO();

    ~GmshIO();

public:
    Mesh read(const eckit::PathName& file_path) const;

    void read(const eckit::PathName& file_path, Mesh& mesh) const;

    /// Write mesh file
    /// Extra file with available mesh information is written to:
    ///  - filename_info.msh
    void write(const Mesh& mesh, const eckit::PathName& file_path) const;

    /// Write field to file
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write(const Field& field, const eckit::PathName& file_path, openmode mode = std::ios::out) const;

    /// Write fieldset to file using FunctionSpace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write(const FieldSet& fieldset, const FunctionSpace&, const eckit::PathName& file_path,
               openmode mode = std::ios::out) const;

    /// Write field to file using Functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write(const Field& field, const FunctionSpace&, const eckit::PathName& file_path,
               openmode mode = std::ios::out) const;

private:
    /// Write fieldset to file using Nodes functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const FieldSet& fieldset, const functionspace::NodeColumns&, const eckit::PathName& file_path,
                        openmode mode = std::ios::out) const;

    /// Write fieldset to file using Nodes functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const FieldSet& fieldset, const functionspace::NoFunctionSpace&,
                        const eckit::PathName& file_path, openmode mode = std::ios::out) const;

    /// Write fieldset to file using Cells functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const FieldSet& fieldset, const functionspace::CellColumns&, const eckit::PathName& file_path,
                        openmode mode = std::ios::out) const;

    /// Write fieldset to file using StructuredColumns functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const FieldSet& fieldset, const functionspace::StructuredColumns&,
                        const eckit::PathName& file_path, openmode mode = std::ios::out) const;

    /// Write field to file using Nodes functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const Field& field, const functionspace::NodeColumns&, const eckit::PathName& file_path,
                        openmode mode = std::ios::out) const;

    /// Write field to file using Cells functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const Field& field, const functionspace::CellColumns&, const eckit::PathName& file_path,
                        openmode mode = std::ios::out) const;

    /// Write field to file using Nodes functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const Field& field, const functionspace::NoFunctionSpace&, const eckit::PathName& file_path,
                        openmode mode = std::ios::out) const;

    /// Write field to file using StructuredColumns functionspace
    ///  Depending on argument "mode", the fields will be appended,
    ///  or existing file will be overwritten
    void write_delegate(const Field& field, const functionspace::StructuredColumns&, const eckit::PathName& file_path,
                        openmode mode = std::ios::out) const;
};

//----------------------------------------------------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

extern "C" {
GmshIO* atlas__Gmsh__new();
void atlas__Gmsh__delete(GmshIO* This);
Mesh::Implementation* atlas__Gmsh__read(GmshIO* This, char* file_path);
void atlas__Gmsh__write(GmshIO* This, Mesh::Implementation* mesh, char* file_path);
Mesh::Implementation* atlas__read_gmsh(char* file_path);
void atlas__write_gmsh_mesh(const Mesh::Implementation* mesh, char* file_path);
void atlas__write_gmsh_fieldset(const field::FieldSetImpl* fieldset, functionspace::FunctionSpaceImpl* function_space,
                                char* file_path, int mode);
void atlas__write_gmsh_field(const field::FieldImpl* field, functionspace::FunctionSpaceImpl* function_space,
                             char* file_path, int mode);
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace output
}  // namespace atlas
