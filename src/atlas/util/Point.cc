/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/util/Point.h"
#include "atlas/util/NormaliseLongitude.h"

#include "atlas/runtime/Exception.h"

namespace atlas {

void PointLonLat::normalise() {
    constexpr util::NormaliseLongitude normalize_from_zero;
    lon() = normalize_from_zero( lon() );
}

void PointLonLat::normalise( double west ) {
    util::NormaliseLongitude normalize_from_west( west );
    lon() = normalize_from_west( lon() );
}

void PointLonLat::normalise( double west, double east ) {
    util::NormaliseLongitude normalize_between_west_and_east( west, east );
    lon() = normalize_between_west_and_east( lon() );
}

// ------------------------------------------------------------------
// C wrapper interfaces to C++ routines

PointXY* atlas__PointXY__new() {
    return new PointXY();
}
void atlas__PointXY__delete( PointXY* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXY" );
    delete This;
}
void atlas__PointXY__print( PointXY* This, std::ostream* channel ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXY" );
    ATLAS_ASSERT( channel != nullptr );
    *channel << *This;
}
PointXY* atlas__PointXY__new_xy(double x, double y) {
    return new PointXY({x, y});
}
double atlas__PointXY__x( PointXY* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXY" );
    return This->x();
}
double atlas__PointXY__y( PointXY* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXY" );
    return This->y();
}
void atlas__PointXY__assign( PointXY* This, double x, double y ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXY" );
    This->assign(x, y);
}

PointXYZ* atlas__PointXYZ__new() {
    return new PointXYZ();
}
void atlas__PointXYZ__delete( PointXYZ* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXYZ" );
    delete This;
}
void atlas__PointXYZ__print( PointXYZ* This, std::ostream* channel ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXYZ" );
    ATLAS_ASSERT( channel != nullptr );
    *channel << *This;
}
PointXYZ* atlas__PointXYZ__new_xyz(double x, double y, double z) {
    return new PointXYZ({x, y, z});
}
double atlas__PointXYZ__x( PointXYZ* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXYZ" );
    return This->x();
}
double atlas__PointXYZ__y( PointXYZ* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXYZ" );
    return This->y();
}
double atlas__PointXYZ__z( PointXYZ* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXYZ" );
    return This->z();
}
void atlas__PointXYZ__assign( PointXYZ* This, double x, double y, double z ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointXYZ" );
    This->assign(x, y, z);
}

PointLonLat* atlas__PointLonLat__new() {
    return new PointLonLat();
}
void atlas__PointLonLat__delete( PointLonLat* This ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    delete This;
}
void atlas__PointLonLat__print( PointLonLat* This, std::ostream* channel ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    ATLAS_ASSERT( channel != nullptr );
    *channel << *This;
}
PointLonLat* atlas__PointLonLat__new_lonlat(double lon, double lat) {
    return new PointLonLat({lon, lat});
}
double atlas__PointLonLat__lon( PointLonLat* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    return This->lon();
}
double atlas__PointLonLat__lat( PointLonLat* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    return This->lat();
}
void atlas__PointLonLat__assign( PointLonLat* This, double lon, double lat ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    This->assign( lon, lat );
}
void atlas__PointLonLat__normalise( PointLonLat* This) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    This->normalise();
}
void atlas__PointLonLat__normalise_west( PointLonLat* This, double west ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    This->normalise( west );
}
void atlas__PointLonLat__normalise_west_east( PointLonLat* This, double west, double east ) {
    ATLAS_ASSERT( This != nullptr, "Cannot access uninitialised atlas_PointLonLat" );
    This->normalise( west, east );
}

// ------------------------------------------------------------------

}  // namespace atlas
