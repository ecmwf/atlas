#ifndef eckit_PointIndex3_h
#define eckit_PointIndex3_h

#include "atlas/Gmsh.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/MeshGen.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Field.hpp"
#include "atlas/Parameters.hpp"

#include "eckit/container/sptree/SPPoint.h"
#include "eckit/container/KDTree.h"
#include "eckit/container/KDMapped.h"
#include "eckit/container/KDMemory.h"

//-----------------------------------------------------------------------------

namespace eckit {

//-----------------------------------------------------------------------------

class KPoint3 : public SPPoint<3> {
public:

    KPoint3(): SPPoint<3>() {}

    double  operator[] (const size_t& i) const { assert(i<3); return x_[i]; }
    double& operator[] (const size_t& i)       { assert(i<3); return x_[i]; }

    KPoint3( double* p ): SPPoint<3>()
    {
        using namespace atlas;
        x_[XX] = p[XX];
        x_[YY] = p[YY];
        x_[ZZ] = p[ZZ];
    }

    KPoint3( const double lat, const double lon ): SPPoint<3>()
    {
        atlas::latlon_to_3d( lat, lon, x_ );
    }

    KPoint3( const double x, const double y, const double z ): SPPoint<3>(x,y)
    {
        x_[atlas::ZZ] = z;
    }

    template < typename T >
    void assign( const T& p )
    {
        using namespace atlas;
        x_[XX] = p[XX];
        x_[YY] = p[YY];
        x_[ZZ] = p[ZZ];
    }

    friend std::ostream& operator<<(std::ostream& s,const KPoint3& p)
    {
        s << '(' << p.x_[atlas::XX] << ","
                 << p.x_[atlas::YY] << ","
                 << p.x_[atlas::ZZ] << ')';
        return s;
    }
};

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
    typedef KPoint3   Point;
    typedef size_t    Payload;
};

//------------------------------------------------------------------------------------------------------

typedef PointKdTree<PointIndex3TreeTrait>  PointIndex3;

//------------------------------------------------------------------------------------------------------

template < typename Tree >
Tree* create_cell_centre_index( atlas::Mesh& mesh )
{
    atlas::FunctionSpace& triags = mesh.function_space( "triags" );
    atlas::FieldT<double>& triags_centres = triags.field<double>( "centre" );

    const size_t npts = triags.bounds()[1];

    std::vector<typename Tree::Value> p;
    p.reserve(npts);

    for( size_t ip = 0; ip < npts; ++ip )
    {
        p.push_back( typename Tree::Value( typename Tree::Point( triags_centres(atlas::XX,ip),
                                                                 triags_centres(atlas::YY,ip),
                                                                 triags_centres(atlas::ZZ,ip) ), ip ) );
    }

    Tree* tree = new Tree();

    tree->build(p.begin(), p.end());

    return tree;
}

//---------------------------------------------------------------------------------------------------------

bool points_equal( const KPoint3& a, const KPoint3& b );

//---------------------------------------------------------------------------------------------------------

} // namespace eckit

#endif

