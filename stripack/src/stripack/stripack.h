/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <vector>
#include <iostream>

namespace stripack {

class Triangulation {
 public:
  Triangulation(size_t N, const double lon[], const double lat[], int lon_stride, int lat_stride, bool reshuffle = true);
  Triangulation(size_t N, const double lon[], const double lat[], bool reshuffle = true);
  Triangulation(size_t N, const double lonlat[], bool reshuffle = true);

  template <typename LonLatVector>
  Triangulation(const LonLatVector& lonlat, bool reshuffle = true):
    Triangulation(lonlat.size(),reinterpret_cast<const double*>(lonlat.data()),reshuffle) {}


  /// Returns true if there is a containing triangle; false if point is outside triangulation
  bool containingTriangleAndBarycentricCoords(
      const std::array<double, 3> & coords,
      int guessIndex,
      std::array<int, 3> & vertexIndices,
      std::array<double, 3> & barycentricCoords) const;

  std::vector<std::array<int,3>> triangles() const;

 public:
  // Triangulation data to use in C++
  size_t num_nodes_;
  std::vector<double> xs_;
  std::vector<double> ys_;
  std::vector<double> zs_;
  std::vector<size_t> random_permutation_;
  std::vector<size_t> inverse_random_permutation_;
  bool reshuffled_;

  // Allocations for STRIPACK connectivity data; should not read from C++
  std::vector<int> list_;
  std::vector<int> lptr_;
  std::vector<int> lend_;
};

}  // namespace stripack
