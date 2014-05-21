/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

// ==================================================================================
// gribs use the following convention: (from Shahram)
//
// Horizontally:  Points scan in the +i (+x) direction
// Vertically:    Points scan in the -j (-y) direction
//
// The way I verified this was to look at our SAMPLE files (which IFS uses).
// I also verified that IFS does not modify the scanning modes
// so whatever the samples say, is the convention
// ==================================================================================

#include "eckit/log/Log.h"
#include "atlas/grid/ReducedGaussianGrid.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------
// Area: Do we check the area.
// Area: Can we assume area is multiple of the grids ? No.
// NPoints: is this just the grids points, or includes area points, of area does not fit grid
//          assumes it is grid points inclusive of the area.

ReducedGaussianGrid::ReducedGaussianGrid()
: gaussianNumber_(0)
{
}

ReducedGaussianGrid::ReducedGaussianGrid( const std::string& hash,
                     const BoundBox& bbox,
                     const std::vector< Point >& pts,
                     const std::vector<double>& latitudes,
                     long gaussianNumber)
: hash_(hash),
  bbox_(bbox),
  points_(pts),
  latitudes_(latitudes),
  gaussianNumber_(gaussianNumber)
{
   Log::info() << "ReducedGaussianGrid" << std::endl;
}

ReducedGaussianGrid::~ReducedGaussianGrid()
{
    Log::info() << "Destroy a ReducedGaussianGrid" << std::endl;
}

void ReducedGaussianGrid::coordinates( Grid::Coords& r ) const
{
    ASSERT( r.size() == points_.size() );

    for( size_t i = 0; i < points_.size(); ++i )
    {
        r.lat(i) = points_[i].lat();
        r.lon(i) = points_[i].lon();
    }
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
