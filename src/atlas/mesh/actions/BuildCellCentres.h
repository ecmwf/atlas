/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_actions_BuildCellCentres_h
#define atlas_actions_BuildCellCentres_h

namespace atlas {

class Mesh;

namespace mesh {
namespace actions {

/// Generates the cell centres on each cell
class BuildCellCentres {
public:

    /// @note Correct only for Linear Triangles and Quadrilaterals
    void operator()(Mesh&) const;

};

} // namespace actions
} // namespace mesh
} // namespace atlas

#endif
