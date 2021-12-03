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

#include "atlas/array/ArrayView.h"
#include "atlas/library/config.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/MicroDeg.h"
#include "atlas/util/Point.h"

namespace atlas {
namespace util {

// ------------------------------------------------------------------------------------

// Forward declarations
class LonLatMicroDeg;
uidx_t unique_lonlat(const LonLatMicroDeg&);

// ------------------------------------------------------------------------------------

/// @brief LonLatMicroDegrees
///
/// Data is stored and accessed in microdegrees in int type (32 bits per
/// coordinate component)
/// This structure is used to compare 2 lon-lat points in sorting algorithms
/// or to use in sending in MPI buffers using minimal memory,
/// The maximum range in degrees is [-2147.48 2147.48]
class LonLatMicroDeg {
public:
    // -- Constructors taking already microdegrees
    LonLatMicroDeg(int lonlat[2]) {
        p_[LON] = lonlat[LON];
        p_[LAT] = lonlat[LAT];
    }
    LonLatMicroDeg(int lon, int lat) {
        p_[LON] = lon;
        p_[LAT] = lat;
    }

    // -- Constructors taking degrees
    LonLatMicroDeg(const double& lon, const double& lat) {
        p_[LON] = microdeg(lon);
        p_[LAT] = microdeg(lat);
    }
    LonLatMicroDeg(const double lonlat[2]) {
        p_[LON] = microdeg(lonlat[LON]);
        p_[LAT] = microdeg(lonlat[LAT]);
    }
    // LonLatMicroDeg( const array::ArrayView<double,1>& lonlat )     {
    // p_[LON]=microdeg(lonlat[LON]); p_[LAT]=microdeg(lonlat[LAT]); }
    LonLatMicroDeg(const PointLonLat& p) {
        p_[LON] = microdeg(p.lon());
        p_[LAT] = microdeg(p.lat());
    }

    // -- Methods
    int lon() const { return p_[LON]; }
    int lat() const { return p_[LAT]; }
    void set_lon(int lon) { p_[LON] = lon; }
    void set_lat(int lat) { p_[LAT] = lat; }
    const int* data() const { return p_; }
    uidx_t unique() { return unique_lonlat(*this); }

    // -- Comparison operator
    bool operator<(const LonLatMicroDeg& other) const;

private:
    int p_[2];
};

// ------------------------------------------------------------------------------------

inline bool LonLatMicroDeg::operator<(const LonLatMicroDeg& other) const {
    if (p_[LAT] > other.p_[LAT])
        return true;
    if (p_[LAT] == other.p_[LAT])
        return (p_[LON] < other.p_[LON]);
    return false;
}

// ------------------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas

#include "atlas/util/Unique.h"
