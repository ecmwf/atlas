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


class CGALSphericalTriangulation {
public:
    template <typename LonLatVector>
    CGALSphericalTriangulation(const LonLatVector& lonlat);
    CGALSphericalTriangulation(size_t N, const double lonlat[]);
    CGALSphericalTriangulation(size_t N, const double lon[], const double lat[]);
    CGALSphericalTriangulation(size_t N, const double lon[], const double lat[], int lon_stride, int lat_stride);

    ~CGALSphericalTriangulation();

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
    struct CGAL;
    std::unique_ptr<CGAL> cgal_;
    std::vector<std::array<double,3>> points_xyz_;
};

template <typename LonLatVector>
inline CGALSphericalTriangulation::CGALSphericalTriangulation(const LonLatVector& lonlat) :
        CGALSphericalTriangulation(lonlat.size(), reinterpret_cast<const double*>(lonlat.data())) {}

template<typename Value>
inline std::vector<std::array<Value,3>> CGALSphericalTriangulation::triangles() const {
    std::vector<std::array<Value,3>> vector(size());
    triangles(vector.data());
    return vector;
}


}
}

