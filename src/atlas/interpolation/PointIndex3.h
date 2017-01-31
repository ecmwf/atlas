/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_interpolation_PointIndex3_h
#define atlas_interpolation_PointIndex3_h

#include "eckit/container/KDMapped.h"
#include "eckit/container/KDMemory.h"
#include "eckit/container/KDTree.h"
#include "eckit/geometry/Point3.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Parameters.h"
#include "atlas/mesh/Mesh.h"


namespace atlas {
namespace interpolation {

//----------------------------------------------------------------------------------------------------------------------

template<class Traits>
class PointKdTree : public eckit::KDTreeMemory<Traits> {
public:
    typedef eckit::KDTreeMemory<Traits> Tree;

    typedef typename Tree::Alloc       Alloc;
    typedef typename Tree::NodeInfo    NodeInfo;
    typedef typename Tree::NodeList    NodeList;
    typedef typename Tree::PayloadType Payload;
    typedef typename Tree::Point       Point;
    typedef typename Tree::Value       Value;

    PointKdTree() : eckit::KDTreeMemory<Traits>() {}
    PointKdTree(Alloc& alloc) : Tree(alloc) {}

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
    typedef size_t                  Payload;
};

typedef PointKdTree<PointIndex3TreeTrait> PointIndex3;

//----------------------------------------------------------------------------------------------------------------------

struct ElemIndex3TreeTrait {
    typedef eckit::geometry::Point3 Point;
    typedef size_t                  Payload;
};

typedef PointKdTree<ElemIndex3TreeTrait> ElemIndex3;

ElemIndex3* create_element_centre_index(const mesh::Mesh& mesh);

//----------------------------------------------------------------------------------------------------------------------

} // namespace interpolation
} // namespace atlas


#endif
