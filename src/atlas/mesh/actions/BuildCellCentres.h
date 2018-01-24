/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

namespace atlas {

class Mesh;
class Field;

namespace mesh {
namespace actions {

/// Generates the cell centres on each cell
class BuildCellCentres {
public:

    BuildCellCentres( const std::string& field_name = "centre" );

    /// @note Correct only for Linear Triangles and Quadrilaterals
    Field& operator()(Mesh&) const;

private:
    std::string field_name_;
};

} // namespace actions
} // namespace mesh
} // namespace atlas
