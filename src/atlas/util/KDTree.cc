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
void atlas__IndexKDTree__reserve( IndexKDTree* This, idx_t size ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return This->reserve(size);
}
void atlas__IndexKDTree__insert( IndexKDTree* This, const Point3* p, const idx_t index ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return This->insert(*p, index);
}
void atlas__IndexKDTree__build( IndexKDTree* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return This->build();
}
void atlas__IndexKDTree__closestPoints( IndexKDTree* This, const Point3* p, size_t k,
                                        Point3*& points, idx_t*& indices, double*& distances) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    ValueList vl(This->closestPoints(*p, k));
    std::vector<Point3> pts;
    std::vector<idx_t> idx;
    std::vector<double> dist;
    for ( size_t i = 0; i < k; ++i ) {
       pts.push_back(vl[i].point());
       idx.push_back(vl[i].payload());
       dist.push_back(vl[i].distance());
    }
    points = new Point3[k];
    indices = new idx_t[k];
    distances = new double[k];
    for ( size_t i = 0; i < k; ++i ) {
       points[i] = pts[i];
       indices[i] = idx[i];
       distances[i] = dist[i];
    }
}
void atlas__IndexKDTree__closestPoint( IndexKDTree* This, const Point3* p,
                                       Point3*& point, idx_t& index, double& distance ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    Value v(This->closestPoint(*p));
    point = new Point3(v.point());
    index = v.payload();
    distance = v.distance();
}
void atlas__IndexKDTree__closestPointsWithinRadius( IndexKDTree* This, const Point3* p, double radius,
                                                    Point3*& points, idx_t*& indices, double*& distances,
                                                    size_t& k) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    ValueList vl(This->closestPoints(*p, radius));
    k = vl.size();
    std::vector<Point3> pts;
    std::vector<idx_t> idx;
    std::vector<double> dist;
    for ( size_t i = 0; i < k; ++i ) {
       pts.push_back(vl[i].point());
       idx.push_back(vl[i].payload());
       dist.push_back(vl[i].distance());
    }
    points = new Point3[k];
    indices = new idx_t[k];
    distances = new double[k];
    for ( size_t i = 0; i < k; ++i ) {
       points[i] = pts[i];
       indices[i] = idx[i];
       distances[i] = dist[i];
    }
}
Geometry* atlas__IndexKDTree__geometry( IndexKDTree* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_IndexKDTree" );
    return new Geometry(This->geometry());
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
