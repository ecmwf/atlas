/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>

#include "eckit/log/Log.h"
#include "atlas/grid/GridSpec.h"

#include "atlas/grid/Gaussian.h"
#include "atlas/grid/Latitudes.h"

using namespace eckit;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

Gaussian::Gaussian( size_t resolution, const BoundBox& bb) :
    resolution_(resolution),
    bound_box_(bb)
{
    ASSERT( resolution_ > 0 );

//    Log::info() << "Build a Gaussian with resolution " << resolution_ << std::endl;

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
    double dlon = ( bb.top_right_.lon() - bb.bottom_left_.lon() ) / nlon ;

    // fill out a std::vector
    std::vector<double> lons;
    double plon = bb.bottom_left_.lon();
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
           coordinates_.push_back( Point(  lats[i], lons[j] ));
           coordinates_.push_back( Point( -lats[i], lons[j] ));
       }
    }
}

Gaussian::~Gaussian()
{
}

std::string Gaussian::hash() const
{
    NOTIMP;
}

Grid::BoundBox Gaussian::boundingBox() const
{
    return bound_box_;
}

void Gaussian::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == coordinates_.size() );

    for( size_t i = 0; i < coordinates_.size(); ++i )
    {
        r.lat(i) = coordinates_[i].lat();
        r.lon(i) = coordinates_[i].lon();
    }
}

GridSpec* Gaussian::spec() const
{
   return new GridSpec(gridType(),"GG");
}


//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
