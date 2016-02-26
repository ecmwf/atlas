/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_PointIndex3_h
#define atlas_PointIndex3_h

#include "eckit/container/KDTree.h"
#include "eckit/container/KDMapped.h"
#include "eckit/container/KDMemory.h"
#include "eckit/geometry/Point3.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/internals/Parameters.h"
#include "atlas/util/io/Gmsh.h"

namespace atlas {
namespace internals {

//----------------------------------------------------------------------------------------------------------------------

template<class Traits>
class PointKdTree : public eckit::KDTreeMemory<Traits> {
public:
    typedef eckit::KDTreeMemory<Traits> Tree;
    typedef typename Tree::Point Point;
    typedef typename Tree::Alloc    Alloc;

    typedef typename Tree::NodeList    NodeList;
    typedef typename Tree::NodeInfo    NodeInfo;
    typedef typename Tree::PayloadType Payload;
    typedef typename Tree::Value   Value;

    PointKdTree(): eckit::KDTreeMemory<Traits>() {}

    NodeInfo nearestNeighbour(const Point& p)
    {
        return Tree::nearestNeighbour(p);
    }

    NodeInfo nearestNeighbourBruteForce(const Point& p)
    {
        return Tree::nearestNeighbourBruteForce(p);
    }

    NodeList findInSphere(const Point& p,double radius)
    {
        return Tree::findInSphere(p, radius);
    }

    NodeList findInSphereBruteForce(const Point& p,double radius)
    {
        return Tree::findInSphereBruteForce(p,radius);
    }

    NodeList kNearestNeighbours(const Point& p,size_t k)
    {
        return Tree::kNearestNeighbours(p, k);
    }

    NodeList kNearestNeighboursBruteForce(const Point& p,size_t k)
    {
        return Tree::kNearestNeighboursBruteForce(p,k);
    }

    PointKdTree(Alloc& alloc): Tree(alloc) {}
};

//----------------------------------------------------------------------------------------------------------------------

struct PointIndex3TreeTrait
{
    typedef eckit::geometry::Point3 Point;
    typedef size_t                  Payload;
};

typedef PointKdTree<PointIndex3TreeTrait>  PointIndex3;

//----------------------------------------------------------------------------------------------------------------------

struct ElemPayload
{
    enum ElementTypeEnum { UNDEFINED=-1, QUAD=0, TRIAG=1 };
    ElemPayload(size_t id, ElementTypeEnum t) : id_(id), type_(t) {}

    size_t id_;
    ElementTypeEnum type_;

    void print(std::ostream& s) const { s << "ElemPayload[id=" << id_ << "]"; }

    bool triangle()      { return type_ == TRIAG; }
    bool quadrilateral() { return type_ == QUAD;  }

    friend std::ostream& operator<<(std::ostream& s, const ElemPayload& p) {
      p.print(s);
      return s;
    }
};

ElemPayload make_elem_payload(size_t id);

struct ElemIndex3TreeTrait
{
    typedef eckit::geometry::Point3 Point;
    typedef ElemPayload             Payload;
};

typedef PointKdTree<ElemIndex3TreeTrait>  ElemIndex3;

ElemIndex3* create_element_centre_index( const atlas::mesh::Mesh& mesh );

//----------------------------------------------------------------------------------------------------------------------

} // namespace internals
} // namespace atlas

#endif

