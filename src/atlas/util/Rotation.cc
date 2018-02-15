/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Rotation.h"

#include <cmath>
#include <iostream>
#include "atlas/util/Constants.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/Earth.h"
#include "eckit/config/Parametrisation.h"

// Temporary option to activate implementation by RMI during ESCAPE
#define OLD_IMPLEMENTATION 0

namespace atlas {
namespace util {

void Rotation::print( std::ostream& out ) const {
    out << "north_pole:" << npole_ << ", south_pole:" << spole_ << ", rotation_angle:" << angle_;
}
std::ostream& operator<<( std::ostream& out, const Rotation& r ) {
    r.print( out );
    return out;
}

PointLonLat Rotation::rotate( const PointLonLat& p ) const {
    PointLonLat rotated( p );
    rotate( rotated.data() );
    return rotated;
}

PointLonLat Rotation::unrotate( const PointLonLat& p ) const {
    PointLonLat unrotated( p );
    unrotate( unrotated.data() );
    return unrotated;
}

void Rotation::precompute() {
    const double theta = -( 90.0 + spole_.lat() ) * Constants::degreesToRadians();
    const double phi   = -spole_.lon() * Constants::degreesToRadians();

    const double sin_theta = std::sin( theta );
    const double cos_theta = std::cos( theta );
    const double sin_phi   = std::sin( phi );
    const double cos_phi   = std::cos( phi );

    // Pt = Rot(z) * Rot(y) * P,   rotate about y axes then z
    // Since we're undoing the rotation described in the definition
    // of the coordinate system,
    // we first rotate by ϑ = -(90 + spole_.lat()) around the y axis
    // (along the rotated Greenwich meridian)
    // and then by φ = -spole_.lon() degrees around the z axis):
    // (xt)   ( cos(φ), sin(φ), 0) (  cos(ϑ), 0, sin(ϑ)) (x)
    // (yt) = (-sin(φ), cos(φ), 0).(  0     , 1, 0     ).(y)
    // (zt)   ( 0     , 0     , 1) ( -sin(ϑ), 0, cos(ϑ)) (z)

    // Expanded
    // xt =  cos(ϑ) cos(φ) x + sin(φ) y + sin(ϑ) cos(φ) z
    // yt = -cos(ϑ) sin(φ) x + cos(φ) y - sin(ϑ) sin(φ) z
    // zt = -sin(ϑ)        x            + cos(ϑ)        z

    rotate_ = {cos_theta * cos_phi,  sin_phi, sin_theta * cos_phi,
               -cos_theta * sin_phi, cos_phi, -sin_theta * sin_phi,
               -sin_theta,           0.,      cos_theta};

    // Assume right hand rule, rotate about z axes and then y
    // P = Rot(y) * Rot(z) * Pt
    // x   (  cos(ϑ), 0, -sin(ϑ)) ( cos(φ), -sin(φ), 0) (xt)
    // y = (  0     , 1,  0     ) ( sin(φ), cos(φ),  0) (yt)
    // z   ( sin(ϑ), 0,   cos(ϑ)) ( 0     , 0     ,  1) (zt)

    // Expanded
    // x   ( cos(ϑ)cos(φ) , -cos(ϑ)sin(φ) , -sin(ϑ)) (xt)
    // y = ( sin(φ)       ,  cos(φ)       ,  0     ).(yt)
    // z   ( sin(ϑ) cos(φ), -sin(ϑ) sin(φ),  cos(ϑ)) (zt)

    unrotate_ = {cos_theta * cos_phi, -cos_theta * sin_phi, -sin_theta, sin_phi, cos_phi, 0.,
                 sin_theta * cos_phi, -sin_theta * sin_phi, cos_theta};

    if ( spole_.lon() == 0. && spole_.lat() == -90. && angle_ == 0. ) rotated_ = false;
    rotation_angle_only_ = spole_.lon() == 0. && spole_.lat() == -90. && rotated_;

    double latrp = ( 90.0 - npole_.lat() ) * Constants::degreesToRadians();
    cos_latrp_   = std::cos( latrp );
    sin_latrp_   = std::sin( latrp );
}

Rotation::Rotation( const PointLonLat& south_pole, double rotation_angle ) {
    spole_ = south_pole;
    npole_ = PointLonLat( spole_.lon() - 180., spole_.lat() + 180. );
    if ( npole_.lat() > 90 ) npole_.lon() += 180.;
    angle_ = rotation_angle;

    precompute();
}

Rotation::Rotation( const eckit::Parametrisation& p ) {
#if OLD_IMPLEMENTATION
    npole_ = {0., 90.};
#endif

    // get rotation angle
    p.get( "rotation_angle", angle_ );

    // get pole
    std::vector<double> pole( 2 );
    if ( p.get( "north_pole", pole ) ) {
        npole_ = PointLonLat( pole.data() );
        spole_ = PointLonLat( npole_.lon() + 180., npole_.lat() - 180. );
        if ( spole_.lat() < -90 ) spole_.lon() -= 180.;
    }
    else if ( p.get( "south_pole", pole ) ) {
        spole_ = PointLonLat( pole.data() );
        npole_ = PointLonLat( spole_.lon() - 180., spole_.lat() + 180. );
        if ( npole_.lat() > 90 ) npole_.lon() += 180.;
    }

    precompute();
}

void Rotation::rotate_old( double crd[] ) const {
    // coordinates of the point P on a rotated sphere with specified pole
    double lon, lat, lont, latt;
    double xt, yt, zt, x, y, z;
    double cos_lon, sin_lon, cos_lat, sin_lat;

    lon     = crd[LON] * Constants::degreesToRadians();
    lat     = crd[LAT] * Constants::degreesToRadians();
    cos_lon = std::cos( lon );
    cos_lat = std::cos( lat );
    sin_lon = std::sin( lon );
    sin_lat = std::sin( lat );

    // cartesian coordinates
    x = cos_lon * cos_lat;
    y = sin_lon * cos_lat;
    z = sin_lat;

    // tilt
    xt = cos_latrp_ * x + sin_latrp_ * z;
    yt = y;
    zt = -sin_latrp_ * x + cos_latrp_ * z;

    // back to spherical coordinates
    lont = std::atan2( yt, xt ) * Constants::radiansToDegrees();
    latt = std::asin( zt ) * Constants::radiansToDegrees();

    // rotate
    crd[0] = lont + npole_.lon();
    crd[1] = latt;
}

void Rotation::unrotate_old( double crd[] ) const {
    // inverse operation of Projection::rotate

    double lont, latt;
    double xt, yt, zt, x, y, z;
    double cos_lont, sin_lont, cos_latt, sin_latt;

    // unrotate
    lont = ( crd[LON] - npole_.lon() ) * Constants::degreesToRadians();
    latt = ( crd[LAT] ) * Constants::degreesToRadians();

    cos_lont = std::cos( lont );
    cos_latt = std::cos( latt );
    sin_lont = std::sin( lont );
    sin_latt = std::sin( latt );

    // cartesian coordinates
    xt = cos_lont * cos_latt;
    yt = sin_lont * cos_latt;
    zt = sin_latt;

    // untilt
    x = cos_latrp_ * xt - sin_latrp_ * zt;
    y = yt;
    z = sin_latrp_ * xt + cos_latrp_ * zt;

    // back to spherical coordinates
    crd[0] = std::atan2( y, x ) * Constants::radiansToDegrees();
    crd[1] = std::asin( z ) * Constants::radiansToDegrees();
}

using RotationMatrix = std::array<std::array<double, 3>, 3>;

inline PointXYZ rotate_geocentric( const PointXYZ& p, const RotationMatrix& R ) {
    return PointXYZ( R[XX][XX] * p.x() + R[XX][YY] * p.y() + R[XX][ZZ] * p.z(),
                     R[YY][XX] * p.x() + R[YY][YY] * p.y() + R[YY][ZZ] * p.z(),
                     R[ZZ][XX] * p.x() + R[ZZ][YY] * p.y() + R[ZZ][ZZ] * p.z() );
}

void Rotation::rotate( double crd[] ) const {
#if OLD_IMPLEMENTATION
    rotate_old( crd );
    return;
#endif

    if ( !rotated_ ) { return; }
    else if ( rotation_angle_only_ ) {
        crd[LON] -= angle_;
        return;
    }

    const PointLonLat L( crd );
    const PointXYZ P     = Earth::convertGeodeticToGeocentric( L, 1. );
    const PointXYZ Pt    = rotate_geocentric( P, rotate_ );
    const PointLonLat Lt = Earth::convertGeocentricToGeodetic( Pt, 1. );

    crd[LON] = Lt.lon() - angle_;
    crd[LAT] = Lt.lat();
}

void Rotation::unrotate( double crd[] ) const {
#if OLD_IMPLEMENTATION
    unrotate_old( crd );
    return;
#endif

    if ( !rotated_ ) { return; }
    else if ( rotation_angle_only_ ) {
        crd[LON] += angle_;
        return;
    }

    PointLonLat Lt( crd );
    Lt.lon() += angle_;

    const PointXYZ Pt   = Earth::convertGeodeticToGeocentric( Lt, 1. );
    const PointXYZ P    = rotate_geocentric( Pt, unrotate_ );
    const PointLonLat L = Earth::convertGeocentricToGeodetic( P, 1. );

    crd[LON] = L.lon();
    crd[LAT] = L.lat();
}

}  // namespace util
}  // namespace atlas
