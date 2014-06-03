/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "PointSet.h"

//------------------------------------------------------------------------------------------------------

#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Parameters.hpp"

using namespace eckit;

namespace atlas {

//------------------------------------------------------------------------------------------------------

PointSet::PointSet( const std::vector< Point >& ipts ) : npts_(ipts.size())
{
    build(ipts);
}

PointSet::PointSet( atlas::Mesh& mesh )
{
    ASSERT( mesh.has_function_space("nodes") );

    FunctionSpace& nodes = mesh.function_space( "nodes" );

    npts_ = nodes.extents()[0];

    ASSERT( npts_ > 0 );

    ASSERT( nodes.has_field("coordinates") );

    ArrayView<double,2> coords ( nodes.field("coordinates") );

    std::vector< PointIndex3::Value > pidx;
    pidx.reserve(npts_);

    for( size_t ip = 0; ip < npts_; ++ip )
        pidx.push_back( PointIndex3::Value( PointIndex3::Point( coords[ip].data() ) , ip ) );

    tree_ = new PointIndex3();

    tree_->build(pidx.begin(), pidx.end());
}

size_t PointSet::search_unique( const Point& p, size_t idx, u_int32_t n  )
{
    PointIndex3::NodeList nearest = tree_->kNearestNeighbours( p, Kn(n) );

    std::vector<size_t> equals;
    equals.reserve( nearest.size() );

    for( size_t i = 0; i < nearest.size(); ++i )
    {
        Point np  = nearest[i].value().point();
        size_t nidx = nearest[i].value().payload();

//            std::cout << "      - " << nidx << " " << np << std::endl;

        if( eckit::geometry::points_equal(p,np) )
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

} // namespace atlas

