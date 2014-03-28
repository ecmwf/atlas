#ifndef atlas_grid_PointIndex3_h
#define atlas_grid_PointIndex3_h

#include "eckit/container/KDTree.h"
#include "eckit/container/KDMapped.h"
#include "eckit/container/KDMemory.h"

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Field.hpp"
#include "atlas/Parameters.hpp"

#include "atlas/grid/Point3.h"
#include "atlas/grid/TriangleIntersection.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {

//------------------------------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------------------------------

struct PointIndex3TreeTrait
{
    typedef Point3   Point;
    typedef size_t    Payload;
};

//------------------------------------------------------------------------------------------------------

typedef PointKdTree<PointIndex3TreeTrait>  PointIndex3;

//------------------------------------------------------------------------------------------------------

PointIndex3* create_cell_centre_index( atlas::Mesh& mesh );

//---------------------------------------------------------------------------------------------------------

} // namespace atlas

#endif

