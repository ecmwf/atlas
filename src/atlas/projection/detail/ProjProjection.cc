/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include "ProjProjection.h"

#include <proj.h>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/Hash.h"

#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Point.h"


namespace atlas {
namespace projection {
namespace detail {


ProjProjection::ProjProjection( const eckit::Parametrisation& param ) :
    sourceToTarget_( nullptr ),
    sourceToGeocentric_( nullptr ),
    context_( PJ_DEFAULT_CTX ) {
    ATLAS_ASSERT( param.get( "proj", proj_ ) && !proj_.empty() );
    param.get( "proj_source", source_ = "EPSG:4326" );          // WGS 84 (lat, lon)
    param.get( "proj_geocentric", geocentric_ = "EPSG:4978" );  // WGS 84 (x, y, z)


    // set x/y transformations to/from lon/lat and to/from geocentric coordinates
    pj_t p1( proj_create_crs_to_crs( context_, source_.c_str(), proj_.c_str(), nullptr ) );
    ATLAS_ASSERT( p1 );
    sourceToTarget_.reset( proj_normalize_for_visualization( context_, p1 ) );
    ATLAS_ASSERT( sourceToTarget_ );

    if ( !geocentric_.empty() ) {
        pj_t p2( proj_create_crs_to_crs( context_, source_.c_str(), geocentric_.c_str(), nullptr ) );
        ATLAS_ASSERT( p2 );
        sourceToGeocentric_.reset( proj_normalize_for_visualization( context_, p2 ) );
        ATLAS_ASSERT( sourceToGeocentric_ );
    }


    // set semi-major/minor axis
    pj_t source( proj_get_source_crs( context_, sourceToTarget_ ) );
    ATLAS_ASSERT( source );

    pj_t ellps( proj_get_ellipsoid( context_, source ) );
    if ( ellps ) {
        double a;
        double b;
        ATLAS_ASSERT( proj_ellipsoid_get_parameters( context_, ellps, &a, &b, nullptr, nullptr ) );
        ATLAS_ASSERT( 0 < b && b <= a );

        eckit::types::is_approximately_equal( a, b )
            ? extraSpec_.set( "radius", a )
            : extraSpec_.set( "semi_major_axis", a ).set( "semi_minor_axis", b );
    }


    // set units
    pj_t target( proj_get_target_crs( context_, sourceToTarget_ ) );
    ATLAS_ASSERT( target );

    pj_t coord( proj_crs_get_coordinate_system( context_, target ) );
    ATLAS_ASSERT( coord );
    ATLAS_ASSERT( proj_cs_get_axis_count( context_, coord ) > 0 );

    const char* units_c_str;
    if ( proj_cs_get_axis_info( nullptr, coord, 0, nullptr, nullptr, nullptr, nullptr, &units_c_str, nullptr,
                                nullptr ) ) {
        std::string units( units_c_str );
        if ( !units.empty() ) {
            extraSpec_.set( "units", units );
        }
    }
}


void ProjProjection::xy2lonlat( double crd[] ) const {
    PJ_COORD P = proj_coord( crd[XX], crd[YY], 0, 0 );
    P          = proj_trans( sourceToTarget_, PJ_INV, P );

    //    std::memcpy(crd, &P, 2 * sizeof(double));
    crd[LON] = P.enu.e;
    crd[LAT] = P.enu.n;
}


void ProjProjection::lonlat2xy( double crd[] ) const {
    PJ_COORD P = proj_coord( crd[LON], crd[LAT], 0, 0 );
    P          = proj_trans( sourceToTarget_, PJ_FWD, P );

    //    std::memcpy(crd, &P, 2 * sizeof(double));
    crd[XX] = P.xy.x;
    crd[YY] = P.xy.y;
}


PointXYZ ProjProjection::xyz( const PointLonLat& lonlat ) const {
    if ( sourceToGeocentric_ ) {
        PJ_COORD P = proj_coord( lonlat.lon(), lonlat.lat(), 0, 0 );
        P          = proj_trans( sourceToGeocentric_, PJ_FWD, P );
        return {P.xyz.x, P.xyz.y, P.xyz.z};
    }
    return ProjectionImpl::xyz( lonlat );
}


RectangularLonLatDomain ProjProjection::lonlatBoundingBox( const Domain& domain ) const {
    return ProjectionImpl::lonlatBoundingBox( domain );
}


ProjProjection::Spec ProjProjection::spec() const {
    Spec spec;
    spec.set( "type", type() ).set( "proj", proj_ ).set( "proj_source", source_ ).set( "proj_geocentric", geocentric_ );
    return spec | extraSpec_;
}


std::string ProjProjection::units() const {
    return extraSpec_.getString( "units", "" );
}


void ProjProjection::hash( eckit::Hash& h ) const {
    h.add( type() );
    h.add( proj_ );
    h.add( source_ );
    h.add( geocentric_ );
}


namespace {
static ProjectionBuilder<ProjProjection> register_1( ProjProjection::static_type() );
}


}  // namespace detail
}  // namespace projection
}  // namespace atlas
