/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Mar 2016

#include "PointsPerLatitude.h"

#include "eckit/memory/ScopedPtr.h"

#include "atlas/grid/detail/pl/classic_gaussian/N.h"

using eckit::Factory;
using eckit::ScopedPtr;

namespace atlas {
namespace grid {
namespace detail {
namespace pl {
namespace classic_gaussian {

//-----------------------------------------------------------------------------

void points_per_latitude_npole_equator (const size_t N, long nlon[]) {
    std::stringstream Nstream;
    Nstream << N;
    std::string Nstr = Nstream.str();
    if( Factory<PointsPerLatitude>::instance().exists(Nstr) ) {
        ScopedPtr<PointsPerLatitude> pl (
            Factory<PointsPerLatitude>::instance().get(Nstr).create() );
        pl->assign(nlon,N);
    } else {
        throw eckit::BadParameter(
            "gaussian::classic::PointsPerLatitude not available for N"+Nstr,
            Here());
    }
}

//-----------------------------------------------------------------------------

void points_per_latitude_npole_spole   (const size_t N, long nlon[]) {
    points_per_latitude_npole_equator(N,nlon);
    size_t end = 2*N-1;
    for(size_t jlat=0; jlat<N; ++jlat) {
        nlon[end--] = nlon[jlat];
    }
}

//-----------------------------------------------------------------------------

}  // namespace classic_gaussian
}  // namespace pl
}  // namespace detail
}  // namespace grid
}  // namespace atlas

