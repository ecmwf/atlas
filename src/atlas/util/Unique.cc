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

uidx_t unique_lonlat(const double& lon, const double& lat, const PeriodicTransform& transform) {
    std::array<double, 2> lonlat{lon, lat};
    transform(lonlat);
    return unique_lonlat(lonlat);
}


uidx_t UniqueLonLat::operator()(const mesh::Connectivity::Row& elem_nodes, const PeriodicTransform& transform) const {
    std::array<double, 2> centroid;
    centroid[LON] = 0.;
    centroid[LAT] = 0.;
    size_t npts   = elem_nodes.size();
    for (size_t jnode = 0; jnode < npts; ++jnode) {
        centroid[LON] += lonlat(elem_nodes(jnode), LON);
        centroid[LAT] += lonlat(elem_nodes(jnode), LAT);
    }
    centroid[LON] /= static_cast<double>(npts);
    centroid[LAT] /= static_cast<double>(npts);
    transform(centroid);
    return unique_lonlat(centroid);
}

uidx_t atlas::util::UniqueLonLat::operator()(int node, const PeriodicTransform& transform) const {
    std::array<double, 2> _lonlat{lonlat(node, LON), lonlat(node, LAT)};
    transform(_lonlat);
    return unique_lonlat(_lonlat[LON], _lonlat[LAT]);
}


}  // namespace util
}  // namespace atlas
