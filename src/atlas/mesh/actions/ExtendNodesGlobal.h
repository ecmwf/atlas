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

#include <string>

namespace atlas {
class Mesh;
class Grid;
}  // namespace atlas

namespace atlas {
namespace mesh {
namespace actions {

/// Adds virtual nodes to the mesh that aren't contained in the Grid Domain
class ExtendNodesGlobal {
public:
    ExtendNodesGlobal(const std::string& gridname = "O16");
    void operator()(const Grid&, Mesh&) const;

private:
    std::string gridname_;
};

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
