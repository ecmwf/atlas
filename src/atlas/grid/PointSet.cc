#include "PointSet.h"

//------------------------------------------------------------------------------------------------------

#include "atlas/Field.hpp"
#include "atlas/FunctionSpace.hpp"
#include "atlas/Mesh.hpp"
#include "atlas/Parameters.hpp"

using namespace atlas;

namespace eckit {

//------------------------------------------------------------------------------------------------------

PointSet::PointSet( const std::vector< KPoint3 >& ipts ) : npts_(ipts.size())
{
    build(ipts);
}

PointSet::PointSet( atlas::Mesh& mesh )
{
    ASSERT( mesh.has_function_space("nodes") );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    npts_ = nodes.bounds()[1];

    ASSERT( npts_ > 0 );

    ASSERT( nodes.has_field("coordinates") );

    FieldT<double>& coords  = nodes.field<double>("coordinates");

    std::vector< PointIndex3::Value > pidx;
    pidx.reserve(npts_);

    for( size_t ip = 0; ip < npts_; ++ip )
        pidx.push_back( PointIndex3::Value( PointIndex3::Point( coords.slice(ip) ) , ip ) );

    tree_ = new PointIndex3();

    tree_->build(pidx.begin(), pidx.end());
}

size_t PointSet::search_unique( const eckit::KPoint3& p, size_t idx, u_int32_t n  )
{
    PointIndex3::NodeList nearest = tree_->kNearestNeighbours( p, Kn(n) );

    std::vector<size_t> equals;
    equals.reserve( nearest.size() );

    for( size_t i = 0; i < nearest.size(); ++i )
    {
        KPoint3 np  = nearest[i].value().point();
        size_t nidx = nearest[i].value().payload();

//            std::cout << "      - " << nidx << " " << np << std::endl;

        if( points_equal(p,np) )
        {
//                std::cout << "      EQUAL !!" << std::endl;
            equals.push_back(nidx);
      }
        else
            break;
    }

    if( equals.size() == nearest.size() ) /* need to increase the search to find all duplicates of this point */
    {
        return this->search_unique(p,idx,++n);
    }
    else /* stop recursion */
    {
        size_t ret = idx; /* if nothing found return same idx */

        if( equals.size() >= 1 ) /* if an equal point was found return the first one */
            ret = equals[0];

        for( size_t n = 1; n < equals.size(); ++n )
            duplicates_[ equals[n] ] = ret;

        return ret;
    }
}

//------------------------------------------------------------------------------------------------------

} // namespace eckit

