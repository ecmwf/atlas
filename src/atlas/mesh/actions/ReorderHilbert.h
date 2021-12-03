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

#include "atlas/mesh/actions/Reorder.h"

namespace atlas {
namespace mesh {
namespace actions {

//----------------------------------------------------------------------------------------------------------------------

/// Reorder implementation that reorders nodes of a mesh following a Hilbert Space-filling curve.
/// Cells and edges are reordered to follow lowest node index.
///
/// Usage:
///     auto reorder = Reorder{ option::type("hilbert") | config };
///     reorder( mesh );
///
/// The optional extra config can contain:
///
///     - "recursion"    : <int>  (default=30)   // Recursion of hilbert space-filling curve,
///                                             // needs to be large enough to provide unique node indices.
///
///    - "ghost_at_end" : <bool> (default=true) // Determines if ghost nodes should be reordered in between
///                                             // internal nodes or added/remain at the end
class ReorderHilbert : public ReorderImpl {
public:
    ReorderHilbert(const eckit::Parametrisation& config = util::NoConfig());

    std::vector<idx_t> computeNodesOrder(Mesh&) override;

private:
    idx_t recursion_{30};
    bool ghost_at_end_{true};
};

// ------------------------------------------------------------------

}  // namespace actions
}  // namespace mesh
}  // namespace atlas
