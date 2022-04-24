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

#include <cmath>
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/LocalView.h"
#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/util/CoordinateEnums.h"
#include "atlas/util/LonLatMicroDeg.h"
#include "atlas/util/MicroDeg.h"
#include "atlas/util/PeriodicTransform.h"

namespace atlas {
namespace util {

class PeriodicTransform;

// ----------------------------------------------------------------------------

/// @brief Compute unique positive index from lon-lat coordinates in
/// microdegrees
/// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
uidx_t unique_lonlat_microdeg(const int lon, const int lat);

/// @brief Compute unique positive index from lon-lat coordinates in
/// microdegrees
/// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
uidx_t unique_lonlat_microdeg(const int lonlat[]);

/// @brief Compute unique positive index from lon-lat coordinates in
/// microdegrees
/// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
uidx_t unique_lonlat(const LonLatMicroDeg&);

/// @brief Compute unique positive index from lon-lat coordinates in degrees
/// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
uidx_t unique_lonlat(const double& lon, const double& lat);
uidx_t unique_lonlat(const double lonlat[]);
uidx_t unique_lonlat(const array::LocalView<double, 1>& lonlat);
uidx_t unique_lonlat(const std::array<double, 2>& lonlat);

/// @brief Compute unique positive index from lon-lat coordinates in degrees
/// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
uidx_t unique_lonlat(const double& lon, const double& lat, const PeriodicTransform&);

/// @brief Compute unique positive index from lon-lat coordinates in degrees.
/// coordinates are stored in order:
/// [ x1, y1,   x2, y2,   ... ,   xn, yn ]
/// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
uidx_t unique_lonlat(const double elem_lonlat[], size_t npts);

/// @brief Compute unique positive index for a element
/// This class is a functor initialised with the nodes
class UniqueLonLat {
public:
    // UniqueLonLat() {}

    /// @brief Constructor, needs nodes functionspace to cache the lonlat field
    UniqueLonLat(const mesh::Nodes&);

    /// @brief Constructor
    UniqueLonLat(const Mesh&);

    /// @brief Compute unique positive index of a node defined by node index.
    /// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
    uidx_t operator()(int node) const;

    uidx_t operator()(int node, const PeriodicTransform& transform) const;

    /// @brief Compute unique positive index of element defined by node indices.
    /// The assumption is that the elements exist in a lon-lat domain and don't
    //  degenerate to a line.
    /// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
    uidx_t operator()(const mesh::Connectivity::Row& elem_nodes) const;

    /// @brief Compute unique positive index of element defined by node indices.
    /// The assumption is that the elements exist in a lon-lat domain and don't
    //  degenerate to a line.
    /// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
    uidx_t operator()(const mesh::Connectivity::Row& elem_nodes, const PeriodicTransform& transform) const;

    /// @brief Compute unique positive index of element defined by node indices.
    /// The assumption is that the elements exist in a lon-lat domain and don't
    //  degenerate to a line.
    /// @return uidx_t Return type depends on ATLAS_BITS_GLOBAL [32/64] bits
    uidx_t operator()(const int elem_nodes[], size_t npts) const;

    /// @brief update the internally cached lonlat view if the field has changed
    void update();

private:
    const mesh::Nodes* nodes;
    array::ArrayView<const double, 2> lonlat;
};

// ----------------------------------------------------------------------------

/* Inline implementation for maximum performance */

namespace detail {

/// lon and lat arguments in microdegrees
template <typename T>
inline T uniqueT(const int lon, const int lat) {
    return uniqueT<T>(lon, lat);
}
template <>
inline int uniqueT<int>(const int lon, const int lat);
template <>
inline long uniqueT<long>(const int lon, const int lat);

/// @brief unique32 computes 32bit positive unique id
/// max precision is 0.02 degrees
/// Numbering follows ECMWF standard grib ordering (from N-S and W-E)
inline int unique32(const int lon_microdeg, const int lat_microdeg) {
    // Truncate microdegrees to order of degrees (16 bits), and add bits together
    // to 32 bit int.
    int iy = static_cast<int>((180000000 - lat_microdeg) * 5e-5);  // (2*microdeg(90.)-lat)*5e-5
    int ix = static_cast<int>((lon_microdeg + 720000000) * 5e-5);  // (lon+2*microdeg(360.))*5e-5
    iy <<= 17;
    int id = iy | ix;
    return id;
}
template <>
inline int uniqueT<int>(const int lon, const int lat) {
    return unique32(lon, lat);
}

/// @brief unique64 computes 64bit positive unique id
/// max precision is 1 microdegree
/// Numbering follows ECMWF standard grib ordering (from N-S and W-E)
inline long unique64(const int lon_microdeg, const int lat_microdeg) {
    // Truncate microdegrees to (32 bits), and add bits together to 64 bit long.
    long iy = 360000000l - long(lat_microdeg);   // (4*microdeg(90.)-lat)
    long ix = long(lon_microdeg) + 1440000000l;  // (lon+4*microdeg(360.))

    iy <<= 31;
    long id = iy | ix;
    return id;
}
template <>
inline long uniqueT<long>(const int lon, const int lat) {
    return unique64(lon, lat);
}
}  // namespace detail

inline uidx_t unique_lonlat_microdeg(const int lon, const int lat) {
    return detail::uniqueT<uidx_t>(lon, lat);
}

inline uidx_t unique_lonlat_microdeg(const int lonlat[]) {
    return detail::uniqueT<uidx_t>(lonlat[LON], lonlat[LAT]);
}

inline uidx_t unique_lonlat(const LonLatMicroDeg& p) {
    return unique_lonlat_microdeg(p.data());
}

inline uidx_t unique_lonlat(const double& lon, const double& lat) {
    return detail::uniqueT<uidx_t>(microdeg(lon), microdeg(lat));
}

inline uidx_t unique_lonlat(const double lonlat[]) {
    return detail::uniqueT<uidx_t>(microdeg(lonlat[LON]), microdeg(lonlat[LAT]));
}

inline uidx_t unique_lonlat(const std::array<double, 2>& lonlat) {
    return detail::uniqueT<uidx_t>(microdeg(lonlat[LON]), microdeg(lonlat[LAT]));
}

inline uidx_t unique_lonlat(const array::LocalView<double, 1>& lonlat) {
    return unique_lonlat(lonlat.data());
}

inline uidx_t unique_lonlat(const double elem_lonlat[], size_t npts) {
    std::array<double, 2> centroid{0., 0.};
    for (size_t jnode = 0; jnode < npts; ++jnode) {
        centroid[LON] += elem_lonlat[jnode * 2 + LON];
        centroid[LAT] += elem_lonlat[jnode * 2 + LAT];
    }
    centroid[LON] /= static_cast<double>(npts);
    centroid[LAT] /= static_cast<double>(npts);

    return unique_lonlat(centroid);
}

inline UniqueLonLat::UniqueLonLat(const Mesh& mesh):
    nodes(&mesh.nodes()), lonlat(array::make_view<double, 2>(nodes->lonlat())) {
    update();
}

inline UniqueLonLat::UniqueLonLat(const mesh::Nodes& _nodes):
    nodes(&_nodes), lonlat(array::make_view<double, 2>(nodes->lonlat())) {
    update();
}

inline uidx_t UniqueLonLat::operator()(int node) const {
    return unique_lonlat(lonlat(node, LON), lonlat(node, LAT));
}

inline uidx_t UniqueLonLat::operator()(const mesh::Connectivity::Row& elem_nodes) const {
    std::array<double, 2> centroid{0., 0.};
    size_t npts = elem_nodes.size();
    for (size_t jnode = 0; jnode < npts; ++jnode) {
        centroid[LON] += lonlat(elem_nodes(jnode), LON);
        centroid[LAT] += lonlat(elem_nodes(jnode), LAT);
    }
    centroid[LON] /= static_cast<double>(npts);
    centroid[LAT] /= static_cast<double>(npts);

    return unique_lonlat(centroid);
}

inline uidx_t UniqueLonLat::operator()(const int elem_nodes[], size_t npts) const {
    std::array<double, 2> centroid{0., 0.};

    for (size_t jnode = 0; jnode < npts; ++jnode) {
        centroid[LON] += lonlat(elem_nodes[jnode], LON);
        centroid[LAT] += lonlat(elem_nodes[jnode], LAT);
    }
    centroid[LON] /= static_cast<double>(npts);
    centroid[LAT] /= static_cast<double>(npts);

    return unique_lonlat(centroid);
}

inline void UniqueLonLat::update() {
    lonlat = array::make_view<double, 2>(nodes->lonlat());
}

// ----------------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
