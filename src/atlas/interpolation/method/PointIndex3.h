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

#include "eckit/container/KDMapped.h"
#include "eckit/container/KDMemory.h"
#include "eckit/container/KDTree.h"
#include "eckit/geometry/Point3.h"

#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

template <class Traits>
class PointKdTree : public eckit::KDTreeMemory<Traits> {
public:
    typedef eckit::KDTreeMemory<Traits> Tree;

    typedef typename Tree::Alloc Alloc;
    typedef typename Tree::NodeInfo NodeInfo;
    typedef typename Tree::NodeList NodeList;
    typedef typename Tree::PayloadType Payload;
    typedef typename Tree::Point Point;
    typedef typename Tree::Value Value;

    PointKdTree(): eckit::KDTreeMemory<Traits>() {}
    PointKdTree(Alloc& alloc): Tree(alloc) {}

    using Tree::findInSphere;
    using Tree::kNearestNeighbours;
    using Tree::nearestNeighbour;

    using Tree::findInSphereBruteForce;
    using Tree::kNearestNeighboursBruteForce;
    using Tree::nearestNeighbourBruteForce;
};

//----------------------------------------------------------------------------------------------------------------------

struct PointIndex3TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef size_t Payload;
};

typedef PointKdTree<PointIndex3TreeTrait> PointIndex3;

//----------------------------------------------------------------------------------------------------------------------

struct ElemIndex3TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef size_t Payload;
};

typedef PointKdTree<ElemIndex3TreeTrait> ElemIndex3;

ElemIndex3* create_element_kdtree(const Mesh& mesh, const Field& field_centres);

// TODO: remove this function, and use "create_element_kdtree(const Field&)"
// instead.
ElemIndex3* create_element_centre_index(const Mesh& mesh);

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
