/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
class Build2DCellCentres {
public:
    Build2DCellCentres(const std::string& field_name = "centre", bool force_recompute = false);
    Build2DCellCentres(eckit::Configuration&);

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
