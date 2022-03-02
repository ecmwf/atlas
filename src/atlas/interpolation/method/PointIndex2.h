/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/container/KDMapped.h"
#include "eckit/container/KDMemory.h"
#include "eckit/container/KDTree.h"
#include "eckit/geometry/Point2.h"

#include "atlas/field/Field.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/util/CoordinateEnums.h"

namespace atlas {
namespace interpolation {
namespace method {

//----------------------------------------------------------------------------------------------------------------------

template <class Traits>
class Point2KdTree : public eckit::KDTreeMemory<Traits> {
public:
    typedef eckit::KDTreeMemory<Traits> Tree;

    typedef typename Tree::Alloc Alloc;
    typedef typename Tree::NodeInfo NodeInfo;
    typedef typename Tree::NodeList NodeList;
    typedef typename Tree::PayloadType Payload;
    typedef typename Tree::Point Point;
    typedef typename Tree::Value Value;

    Point2KdTree(): eckit::KDTreeMemory<Traits>() {}
    Point2KdTree(Alloc& alloc): Tree(alloc) {}

    using Tree::findInSphere;
    using Tree::kNearestNeighbours;
    using Tree::nearestNeighbour;

    using Tree::findInSphereBruteForce;
    using Tree::kNearestNeighboursBruteForce;
    using Tree::nearestNeighbourBruteForce;
};

//----------------------------------------------------------------------------------------------------------------------

struct PointIndex2TreeTrait {
    typedef eckit::geometry::Point2 Point;
    typedef size_t Payload;
};

typedef Point2KdTree<PointIndex2TreeTrait> PointIndex2;

//----------------------------------------------------------------------------------------------------------------------

struct ElemIndex2TreeTrait {
    typedef eckit::geometry::Point2 Point;
    typedef size_t Payload;
};

typedef Point2KdTree<ElemIndex2TreeTrait> ElemIndex2;

ElemIndex2* create_element2D_kdtree(const Mesh& mesh, const Field& field_centres);

//----------------------------------------------------------------------------------------------------------------------

}  // namespace method
}  // namespace interpolation
}  // namespace atlas
