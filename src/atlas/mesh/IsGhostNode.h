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

#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace mesh {

class IsGhostNode {
public:
    IsGhostNode(const mesh::Nodes& nodes):
        flags_(array::make_view<int, 1>(nodes.flags())), ghost_(array::make_view<int, 1>(nodes.ghost())) {}

    bool operator()(idx_t idx) const { return Nodes::Topology::check(flags_(idx), Nodes::Topology::GHOST); }

private:
    array::ArrayView<const int, 1> flags_;
    array::ArrayView<const int, 1> ghost_;
};

}  // namespace mesh
}  // namespace atlas
