/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/geometry/Point3.h"

#include "atlas/grid.h"
#include "atlas/util/KDTree.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/detail/KDTree.h"

namespace atlas {
namespace util {

using Handle         = typename ObjectHandle<detail::KDTreeBase<idx_t, Point3>>::Handle;
using Implementation = typename Handle::Implementation;
using Value          = typename Implementation::Value;
using ValueList      = typename Implementation::ValueList;

// C wrapper interfaces to C++ routines
IndexKDTree* atlas__IndexKDTree__new() {
    return new IndexKDTree();
}
IndexKDTree* atlas__IndexKDTree__new_geometry( const Geometry* geometry ) {
    return new IndexKDTree(*geometry);
}
void atlas__IndexKDTree__delete( IndexKDTree* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    delete This;
}
void atlas__IndexKDTree__reserve( IndexKDTree* This, const idx_t size ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return This->reserve(size);
}
void atlas__IndexKDTree__insert( IndexKDTree* This, const double lon, const double lat, const idx_t index ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return This->insert(PointLonLat{lon, lat}, index);
}
void atlas__IndexKDTree__build( IndexKDTree* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return This->build();
}
void atlas__IndexKDTree__closestPoints( IndexKDTree* This, const double plon, const double plat, const size_t k,
                                        double*& lons, double*& lats, idx_t*& indices, double*& distances ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    ValueList vl(This->closestPoints(PointLonLat{plon, plat}, k));
    lons = new double[k];
    lats = new double[k];
    indices = new idx_t[k];
    distances = new double[k];
    for ( size_t i = 0; i < k; ++i ) {
      PointLonLat lonlat;
      This->geometry().xyz2lonlat(vl[i].point(), lonlat);
      lonlat.normalise();
      lons[i] = lonlat.lon();
      lats[i] = lonlat.lat();
      indices[i] = vl[i].payload();
      distances[i] = vl[i].distance();
    }
}
void atlas__IndexKDTree__closestPoint( IndexKDTree* This, const double plon, const double plat,
                                       double& lon, double& lat, idx_t& index, double& distance ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    Value v(This->closestPoint(PointLonLat{plon, plat}));
    PointLonLat lonlat;
    This->geometry().xyz2lonlat(v.point(), lonlat);
    lonlat.normalise();
    lon = lonlat.lon();
    lat = lonlat.lat();
    index = v.payload();
    distance = v.distance();
}
void atlas__IndexKDTree__closestPointsWithinRadius( IndexKDTree* This, const double plon, const double plat, double radius,
                                                    size_t& k, double*& lons, double*& lats, idx_t*& indices, double*& distances ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    ValueList vl(This->closestPointsWithinRadius(PointLonLat{plon, plat}, radius));
    k = vl.size();
    lons = new double[k];
    lats = new double[k];
    indices = new idx_t[k];
    distances = new double[k];
    for ( size_t i = 0; i < k; ++i ) {
      PointLonLat lonlat;
      This->geometry().xyz2lonlat(vl[i].point(), lonlat);
      lonlat.normalise();
      lons[i] = lonlat.lon();
      lats[i] = lonlat.lat();
      indices[i] = vl[i].payload();
      distances[i] = vl[i].distance();
    }
}
Geometry* atlas__IndexKDTree__geometry( IndexKDTree* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return new Geometry(This->geometry());
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
