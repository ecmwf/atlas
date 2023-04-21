/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <array>
#include <vector>
#include <memory>

#include "atlas/library/config.h"
namespace atlas {
namespace util {

class QhullSphericalTriangulation {
public:
    template <typename LonLatVector>
    QhullSphericalTriangulation(const LonLatVector& lonlat);
    QhullSphericalTriangulation(size_t N, const double lonlat[]);
    QhullSphericalTriangulation(size_t N, const double lon[], const double lat[]);
    QhullSphericalTriangulation(size_t N, const double lon[], const double lat[], int lon_stride, int lat_stride);

    ~QhullSphericalTriangulation();

    /// @return number of triangles
    size_t size() const;

    void triangles(std::array<int,3>[]) const;
    void triangles(std::array<long,3>[]) const;
    void triangles(std::array<long long,3>[]) const;
    void triangles(std::array<size_t,3>[]) const;

    template<typename Value>
    std::vector<std::array<Value,3>> triangles() const;
    std::vector<std::array<idx_t,3>> triangles() const;

private:
    struct Qhull;
    std::unique_ptr<Qhull> qhull_;
    std::vector<std::array<double,3>> points_xyz_;
};

template <typename LonLatVector>
inline QhullSphericalTriangulation::QhullSphericalTriangulation(const LonLatVector& lonlat) :
        QhullSphericalTriangulation(lonlat.size(), reinterpret_cast<const double*>(lonlat.data())) {}

template<typename Value>
inline std::vector<std::array<Value,3>> QhullSphericalTriangulation::triangles() const {
    std::vector<std::array<Value,3>> vector(size());
    triangles(vector.data());
    return vector;
}


}
}

