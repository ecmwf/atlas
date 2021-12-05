/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Mar 2016

#include "PointsPerLatitude.h"

#include <memory>
#include <sstream>

#include "atlas/grid/detail/pl/classic_gaussian/N.h"
#include "atlas/runtime/Exception.h"


namespace atlas {
namespace grid {
namespace detail {
namespace pl {
namespace classic_gaussian {

//-----------------------------------------------------------------------------

template <typename Int>
void points_per_latitude_npole_equator_impl(const size_t N, Int nlon[]) {
    std::stringstream Nstream;
    Nstream << N;
    std::string Nstr = Nstream.str();
    if (PointsPerLatitudeFactory::has(Nstr)) {
        std::unique_ptr<const PointsPerLatitude> pl(PointsPerLatitudeFactory::build(Nstr));
        pl->assign(nlon, N);
    }
    else {
        throw_Exception("gaussian::classic::PointsPerLatitude not available for N" + Nstr, Here());
    }
}

//-----------------------------------------------------------------------------

template <typename Int>
void points_per_latitude_npole_spole_impl(const size_t N, Int nlon[]) {
    points_per_latitude_npole_equator(N, nlon);
    size_t end = 2 * N - 1;
    for (size_t jlat = 0; jlat < N; ++jlat) {
        nlon[end--] = nlon[jlat];
    }
}

//-----------------------------------------------------------------------------

void points_per_latitude_npole_equator(const size_t N, long nlon[]) {
    points_per_latitude_npole_equator_impl(N, nlon);
}
void points_per_latitude_npole_equator(const size_t N, int nlon[]) {
    points_per_latitude_npole_equator_impl(N, nlon);
}

//-----------------------------------------------------------------------------

void points_per_latitude_npole_spole(const size_t N, long nlon[]) {
    points_per_latitude_npole_spole_impl(N, nlon);
}
void points_per_latitude_npole_spole(const size_t N, int nlon[]) {
    points_per_latitude_npole_spole_impl(N, nlon);
}

//-----------------------------------------------------------------------------

}  // namespace classic_gaussian
}  // namespace pl
}  // namespace detail
}  // namespace grid
}  // namespace atlas
