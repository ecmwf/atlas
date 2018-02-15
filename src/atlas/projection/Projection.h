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

#include "eckit/config/Parametrisation.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Point.h"

//---------------------------------------------------------------------------------------------------------------------

// Forward declarations
namespace eckit {
class Parametrisation;
}

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {

//---------------------------------------------------------------------------------------------------------------------

class Projection {
public:
    using Implementation = projection::detail::ProjectionImpl;
    using Spec           = Implementation::Spec;

public:
    Projection();
    Projection( const Projection& );
    Projection( const Implementation* );
    Projection( const eckit::Parametrisation& );

    void xy2lonlat( double crd[] ) const;
    void lonlat2xy( double crd[] ) const;

    PointLonLat lonlat( const PointXY& ) const;
    PointXY xy( const PointLonLat& ) const;

    bool strictlyRegional() const;

    Spec spec() const;

    std::string units() const;

    operator bool() const;

    std::string type() const { return projection_->type(); }

    void hash( eckit::Hash& ) const;

private:
    eckit::SharedPtr<Implementation> projection_;
};

//---------------------------------------------------------------------------------------------------------------------

inline void Projection::xy2lonlat( double crd[] ) const {
    return projection_->xy2lonlat( crd );
}
inline void Projection::lonlat2xy( double crd[] ) const {
    return projection_->lonlat2xy( crd );
}
inline PointLonLat Projection::lonlat( const PointXY& xy ) const {
    return projection_->lonlat( xy );
}
inline PointXY Projection::xy( const PointLonLat& lonlat ) const {
    return projection_->xy( lonlat );
}
inline bool Projection::strictlyRegional() const {
    return projection_->strictlyRegional();
}
inline Projection::Spec Projection::spec() const {
    return projection_->spec();
}
inline std::string Projection::units() const {
    return projection_->units();
}
inline Projection::operator bool() const {
    return projection_->operator bool();
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
