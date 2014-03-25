/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"

#include "eckit/grid/Gaussian.h"
#include "eckit/grid/Latitudes.h"
#include <math.h>

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {

//-----------------------------------------------------------------------------

Gaussian::Gaussian( size_t resolution, const BoundBox2D& bb) :
    resolution_(resolution),
    bound_box_(bb)
{
    ASSERT( resolution_ > 0 );

    Log::info() << "Build a Gaussian with resolution " << resolution_ << std::endl;

    // Regular Gaussian grids have four times resolution_ along each
    // latitude. The number of latitudes in each hemisphere == resolution_
    // NB There is no ZERO meridian.

    const size_t nlon = resolution_ * 4;
    coordinates_.reserve( (resolution_ * 2) * nlon );

    // we generate latitude values according to a formula. The spacing of the 
    // latitudes is not constant with latitude.
    //
    // with number of lats in the N hemisphere equal to the resolution_
    // with these replicated in the S hemisphere
    // NB we don't go all the way to 90 degrees in either case
    // e.g. http://www.ecmwf.int/publications/manuals/libraries/interpolation/n80FIS.html

    std::vector<double> lats;
    Latitudes::gaussian(resolution_, lats);
    ASSERT(lats.size() == resolution_);

    // generate longitudes common to all latitudes for regular gaussian grids
    // 
    // work out the delta between longitudes
    double dlon = ( bb.top_right_.lon_ - bb.bottom_left_.lon_ ) / nlon ;

    // fill out a std::vector
    std::vector<double> lons;
    double plon = bb.bottom_left_.lon_;
    for (size_t j = 0; j <= nlon; ++j)
    {
        lons.push_back(plon);
        plon += dlon;
    }

    for( size_t i = 0; i < lats.size(); ++i )
    {
        // now compute all longitudes along this latitude
        /// @todo longitudes can be precomputed for regular grids
        /// @todo we wrap longitudes to include 0 and 360 at present
       for (size_t j = 0; j < lons.size(); j++)
       {
           // write to both hemispheres
           coordinates_.push_back( Point2D(  lats[i], lons[j] ));
           coordinates_.push_back( Point2D( -lats[i], lons[j] ));
       }
    }

}

Gaussian::~Gaussian()
{
}

BoundBox2D Gaussian::boundingBox() const
{
    return bound_box_;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
