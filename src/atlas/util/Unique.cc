/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */
 
#include "atlas/util/Unique.h"
#include "atlas/util/PeriodicTransform.h"

namespace atlas {
namespace util {

uidx_t unique_lonlat( const double& lon, const double& lat, const PeriodicTransform& transform ) {
    std::array<double,2> lonlat{lon,lat};
    transform( lonlat );
    return detail::uniqueT<uidx_t>( microdeg( lonlat[LON] ), microdeg( lonlat[LAT] ) );
}

uidx_t UniqueLonLat::operator()( const mesh::Connectivity::Row& elem_nodes, const PeriodicTransform& transform) const {
    double centroid[2];
    centroid[XX] = 0.;
    centroid[YY] = 0.;
    size_t npts  = elem_nodes.size();
    for ( size_t jnode = 0; jnode < npts; ++jnode ) {
        centroid[XX] += xy( elem_nodes( jnode ), XX );
        centroid[YY] += xy( elem_nodes( jnode ), YY );
    }
    centroid[XX] /= static_cast<double>( npts );
    centroid[YY] /= static_cast<double>( npts );
    return unique_lonlat( centroid[XX], centroid[YY], transform );
}



}  // namespace util
}  // namespace atlas
