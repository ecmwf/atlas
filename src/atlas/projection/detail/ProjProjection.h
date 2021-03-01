/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include <memory>

#include "atlas/projection/detail/ProjectionImpl.h"
#include "atlas/util/Config.h"


extern "C" {
using PJ = struct PJconsts;
PJ* proj_destroy( PJ* );
using PJ_CONTEXT = struct projCtx_t;
PJ_CONTEXT* proj_context_destroy( PJ_CONTEXT* );
}


namespace atlas {
namespace projection {
namespace detail {


class ProjProjection final : public ProjectionImpl {
public:
    // -- Types
    // None

    // -- Exceptions
    // None

    // Constructors

    ProjProjection( const eckit::Parametrisation& );

    // Destructor
    // None

    // -- Methods

    static std::string static_type() { return "proj"; }

    // -- Overridden methods

    std::string type() const override { return static_type(); }

    void xy2lonlat( double[] ) const override;
    void lonlat2xy( double[] ) const override;

    Jacobian jacobian( const PointLonLat& ) const override;

    PointXYZ xyz( const PointLonLat& ) const override;

    bool strictlyRegional() const override { return false; }
    RectangularLonLatDomain lonlatBoundingBox( const Domain& ) const override;

    Spec spec() const override;
    std::string units() const override;
    void hash( eckit::Hash& ) const override;

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Members
    // None

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None

public:
    // -- Types

    struct pj_t : std::unique_ptr<PJ, decltype( &proj_destroy )> {
        using t = std::unique_ptr<PJ, decltype( &proj_destroy )>;
        explicit pj_t( PJ* ptr ) : t( ptr, &proj_destroy ) {}
        operator PJ*() const { return t::get(); }
    };

    struct ctx_t : std::unique_ptr<PJ_CONTEXT, decltype( &proj_context_destroy )> {
        using t = std::unique_ptr<PJ_CONTEXT, decltype( &proj_context_destroy )>;
        explicit ctx_t( PJ_CONTEXT* ptr ) : t( ptr, &proj_context_destroy ) {}
        operator PJ_CONTEXT*() const { return t::get(); }
    };

    // -- Contructors
    // None

    // -- Destructor
    // None

    // -- Convertors
    // None

    // -- Operators
    // None

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Members
private:
    Normalise normalise_;
    std::string proj_;
    std::string source_;
    std::string geocentric_;
    bool source_encoded_;
    bool geocentric_encoded_;

    pj_t sourceToTarget_;
    pj_t sourceToGeocentric_;
    ctx_t context_;

    Spec extraSpec_;

    // -- Methods
    // None

    // -- Overridden methods
    // None

    // -- Class members
    // None

    // -- Class methods
    // None

    // -- Friends
    // None (how sad)
};


}  // namespace detail
}  // namespace projection
}  // namespace atlas
