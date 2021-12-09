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

#include "atlas/util/Config.h"

namespace atlas {

class Mesh;
class Field;

namespace mesh {
namespace actions {

/// Generates the cell centres on each cell
class BuildCellCentres {
public:
    BuildCellCentres(const std::string& field_name = "centre", bool force_recompute = false);
    BuildCellCentres(eckit::Configuration&);

    /// @note Correct only for Linear Triangles and Quadrilaterals
    Field& operator()(Mesh&) const;

private:
    std::string field_name_;
    bool force_recompute_;
    bool flatten_virtual_elements_;
};

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
