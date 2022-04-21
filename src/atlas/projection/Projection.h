/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/domain/Domain.h"
#include "atlas/library/config.h"
#include "atlas/projection/Jacobian.h"
#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/ObjectHandle.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
class Hash;
namespace geometry {
class Point2;
}
}  // namespace eckit

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {

class PointLonLat;
class PointXY;
using Point2 = eckit::geometry::Point2;

//---------------------------------------------------------------------------------------------------------------------
namespace util {
class Config;
}
namespace projection {
namespace detail {
class ProjectionImpl;
}
}  // namespace projection

class Projection : DOXYGEN_HIDE(public util::ObjectHandle<projection::detail::ProjectionImpl>) {
public:
    using Spec     = util::Config;
    using Jacobian = projection::Jacobian;

public:
    using Handle::Handle;
    Projection();
    Projection(const eckit::Parametrisation&);
    Projection(const std::string& type, const eckit::Parametrisation&);

    operator bool() const;

    void xy2lonlat(double crd[]) const;
    void xy2lonlat(Point2&) const;
    void lonlat2xy(double crd[]) const;
    void lonlat2xy(Point2&) const;

    Jacobian jacobian(const PointLonLat&) const;

    PointLonLat lonlat(const PointXY&) const;
    PointXY xy(const PointLonLat&) const;

    bool strictlyRegional() const;
    RectangularLonLatDomain lonlatBoundingBox(const Domain&) const;

    Spec spec() const;

    std::string units() const;

    std::string type() const;

    void hash(eckit::Hash&) const;
};

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
