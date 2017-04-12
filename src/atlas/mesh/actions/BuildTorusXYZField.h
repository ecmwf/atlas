/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_actions_BuildTorusXYZField_h
#define atlas_actions_BuildTorusXYZField_h

#include <string>

#include "atlas/grid/Domain.h"
#include "atlas/grid/detail/domain/RectangularDomain.h"

namespace atlas {
  class Field;
}

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
class BuildTorusXYZField {

public:

    explicit BuildTorusXYZField(const std::string& name = "xyz");

    Field& operator()(Mesh&, const atlas::grid::Domain& , double r0, double r1, int nx, int ny) const;
    Field& operator()(mesh::Nodes&,const atlas::grid::Domain& , double r0, double r1, int nx, int ny) const;

private:

    std::string name_;

};

//----------------------------------------------------------------------------------------------------------------------

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif
