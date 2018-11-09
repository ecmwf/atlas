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

#include "atlas/projection/detail/ProjectionImpl.h"

namespace atlas {
namespace projection {
namespace detail {

template <typename Rotation>
class LonLatProjectionT : public ProjectionImpl {
private:
    friend class ProjectionImpl;
    LonLatProjectionT() : ProjectionImpl() {}

public:
    // constructor
    LonLatProjectionT( const eckit::Parametrisation& );

    // destructor
    ~LonLatProjectionT() {}

    // class name
    static std::string static_type() { return Rotation::typePrefix() + "lonlat"; }
    virtual std::string type() const override { return static_type(); }

    // projection and inverse projection
    virtual void xy2lonlat( double crd[] ) const override { rotation_.rotate( crd ); }
    virtual void lonlat2xy( double crd[] ) const override { rotation_.unrotate( crd ); }

    virtual bool strictlyRegional() const override { return false; }

    // specification
    virtual Spec spec() const override;

    virtual std::string units() const override { return "degrees"; }

    virtual operator bool() const override { return rotation_.rotated(); }

    virtual void hash( eckit::Hash& ) const override;

private:
    Rotation rotation_;
};

typedef LonLatProjectionT<NotRotated> LonLatProjection;
typedef LonLatProjectionT<Rotated> RotatedLonLatProjection;

// --------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace projection
}  // namespace atlas
