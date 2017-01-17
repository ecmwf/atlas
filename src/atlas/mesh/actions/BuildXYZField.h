/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_actions_BuildXYZField_h
#define atlas_actions_BuildXYZField_h

#include <string>

namespace atlas {
namespace field {
    class Field;
} }

namespace atlas {
namespace mesh {
  class Mesh;
  class Nodes;
} }

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

/// Creates a XYZ field from the (lon,lat) field
class BuildXYZField {
public:

    explicit BuildXYZField(const std::string& name = "xyz");

    field::Field& operator()(Mesh&) const;
    field::Field& operator()(mesh::Nodes&) const;

private:

    std::string name_;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif
