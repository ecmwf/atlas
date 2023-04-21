/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "stripack.h"

#include <algorithm>
#include <array>
#include <vector>
#include <random>
#include <cassert>
#include <iostream>
#include <string>
#include <ctime>

// Interface for Fortran implementations
extern "C" {
  void stripack_trmesh_f90(const int &, const double[], const double[], const double[],
                           int[], int[], int[], int &, int[], int[], double[], int&);
  void stripack_trfind_f90(const int &, const double[],
                           const int &, const double[], const double[], const double[],
                           const int[], const int[], const int[],
                           double[], int[]);
  void stripack_trlist2_f90(const int &, const int[], const int[], const int[],
                            int &, int[], int &);
};  // extern "C"

namespace stripack {

#define ASSERT assert

class NormaliseLongitude {
public:
    // Normalise longitude between (west - eps, east - eps ) with west = 0., east = 360.
    constexpr NormaliseLongitude(): west_(-eps_), east_(360. - eps_) {}

    // Normalise longitude between ( west-eps, east-eps )  with east = west + 360
    constexpr NormaliseLongitude(double west): west_(west - eps_), east_(west + 360. - eps_) {}

    // Normalise longitude between ( west-eps, east+eps )
    constexpr NormaliseLongitude(double west, double east): west_(west - eps_), east_(east + eps_) {}

    constexpr NormaliseLongitude(const NormaliseLongitude& other): west_(other.west_), east_(other.east_) {}

    double operator()(double lon) const {
        while (lon < west_) {
            lon += 360.;
        }
        while (lon > east_) {
            lon -= 360.;
        }
        return lon;
    }

private:
    double west_;
    double east_;

public:
    static constexpr double eps_ = 1.e-11;
};

void lonlat2xyz(double lon, double lat, double& x, double& y, double& z) {
    /*
     * See https://en.wikipedia.org/wiki/Reference_ellipsoid#Coordinates
     * numerical conditioning for both ϕ (poles) and λ (Greenwich/Date Line).
     *
     * cos α = sqrt( 1 - sin^2 α) is better conditioned than explicit cos α, and
     * coupled with λ in [-180°, 180°[ the accuracy of the trigonometric
     * functions is the same (before converting/multiplying its angle argument
     * to radian) and explicitly chosing -180° over 180° for longitude.
     *
     * These three conditionings combined project very accurately to the sphere
     * poles and quadrants.
     */
    constexpr double degrees_to_radians = M_PI / 180.;
    constexpr NormaliseLongitude normalise_longitude(-180.);
    const double lambda_deg = normalise_longitude(lon);
    const double lambda     = degrees_to_radians * lambda_deg;
    const double phi        = degrees_to_radians * lat;

    const double sin_phi    = std::sin(phi);
    const double cos_phi    = std::sqrt(1. - sin_phi * sin_phi);
    const double sin_lambda = std::abs(lambda_deg) < 180. ? std::sin(lambda) : 0.;
    const double cos_lambda = std::abs(lambda_deg) > 90. ? std::cos(lambda) : std::sqrt(1. - sin_lambda * sin_lambda);

    x = cos_phi * cos_lambda;
    y = cos_phi * sin_lambda;
    z = sin_phi;
}

template<class RandomIt>
void shuffle(RandomIt begin, RandomIt end, unsigned int seed =
             static_cast<std::uint32_t>(std::time(nullptr)), bool reset = false) {
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    typedef std::uniform_int_distribution<diff_t> distr_t;
    typedef typename distr_t::param_type param_t;

    static std::mt19937 generator(seed);
    if (reset) generator.seed(seed);

    distr_t distribution;
    diff_t n = end - begin;
    for (diff_t i = n - 1; i > 0; --i) {
        std::swap(begin[i], begin[distribution(generator, param_t(0, i))]);
    }
}

Triangulation::Triangulation(size_t N, const double lonlat[], bool reshuffle) :
  Triangulation(N, lonlat, lonlat+1, 2, 2, reshuffle) {
}

Triangulation::Triangulation(size_t N, const double lon[], const double lat[], bool reshuffle) :
  Triangulation(N, lon, lat, 1, 1, reshuffle) {
}

Triangulation::Triangulation(size_t N, const double lon[], const double lat[], int lon_stride, int lat_stride, bool reshuffle) :
    num_nodes_(N),
    xs_(num_nodes_),
    ys_(num_nodes_),
    zs_(num_nodes_),
    reshuffled_(false),
    list_(6*(num_nodes_-2)),
    lptr_(6*(num_nodes_-2)),
    lend_(num_nodes_)
{
  // The source latlons come from a model grid so may be highly structured. But triangulating by
  // traversing through structured points can lead to a very non-uniform intermediate triangulation
  // with pathological configurations. So we randomize the source coordinates before triangulation.
  //
  // First, construct a random permutation and its inverse:
  if (reshuffle) {
    random_permutation_.resize(num_nodes_);
    inverse_random_permutation_.resize(num_nodes_);

    std::vector<size_t> indices(num_nodes_);
    std::iota(indices.begin(), indices.end(), 0);
    const unsigned int seed = 7452;   // arbitrary, but reproducible
    shuffle(begin(indices), end(indices), seed);
    for (size_t i = 0; i < num_nodes_; ++i) {
      random_permutation_[i] = indices[i];
      inverse_random_permutation_[indices[i]] = i;
    }
    // Then permute source points while converting from latlon coords to xyz coords
    for (size_t i = 0; i < num_nodes_; ++i) {
      const size_t ip = inverse_random_permutation_[i];
      lonlat2xyz(lon[ip*lon_stride], lat[ip*lat_stride] , xs_[i], ys_[i], zs_[i]);
    }
    reshuffled_ = true;
  }
  else {
    for (size_t i = 0; i < num_nodes_; ++i) {
      const size_t ip = inverse_random_permutation_[i];
      lonlat2xyz(lon[ip*lon_stride], lat[ip*lat_stride] , xs_[i], ys_[i], zs_[i]);
    }
  }


  // Allocations for STRIPACK workspaces
  std::vector<int> near(num_nodes_);
  std::vector<int> next(num_nodes_);
  std::vector<double> dist(num_nodes_);
  int lnew;
  int ier;
  stripack_trmesh_f90(
      static_cast<int>(num_nodes_),
      xs_.data(),   ys_.data(),   zs_.data(),
      list_.data(), lptr_.data(), lend_.data(), lnew,
      near.data(), next.data(), dist.data(), ier);
  if (ier != 0) {
    if (ier==-1) {
      throw std::runtime_error("stripack: Could not triangulate points. Reason: N<3 on input");
    }
    else if (ier==-2) {
      throw std::runtime_error("stripack: Could not triangulate points. Reason: first three nodes are collinear");
    }
    else {
      int l = ier;
      // nodes L and M coincide for some L < M.  The data structure 
      // represents a triangulation of nodes 1 to M-1 in this case.
      if (reshuffled_) {
        l = inverse_random_permutation_[l];
      }
      throw std::runtime_error("stripack: Could not triangulate points. Reason: Node "+std::to_string(l)+
        " is a duplicate.");
    }
    throw std::runtime_error("stripack: Could not triangulate points");
  }
} 

std::vector<std::array<int,3>> Triangulation::triangles() const {
    int ier;
    int num_triag;
    std::vector<std::array<int,3>> triangles(2*num_nodes_-4);
    int* ltri = reinterpret_cast<int*>( triangles.data() );
    stripack_trlist2_f90(
        static_cast<int>(num_nodes_),
        list_.data(), lptr_.data(), lend_.data(), num_triag, ltri, ier);
    if (ier != 0) {
      throw std::runtime_error("stripack: Could infer triangles");
    }

    triangles.resize(num_triag);

    if (reshuffled_) {
        for( size_t j=0; j<3*num_triag; ++j) {
            int& n = ltri[j];
            n = inverse_random_permutation_[n];
        }
    }
    return triangles;
}


bool Triangulation::containingTriangleAndBarycentricCoords(
      const std::array<double, 3> & coords,
      const int guess_index,
      std::array<int, 3> & vertexIndices,
      std::array<double, 3> & barycentricCoords) const {
    ASSERT(guess_index >= 0 && guess_index <= static_cast<int>(num_nodes_) - 1);

    int idx = guess_index;
    if (reshuffled_) {
      idx = random_permutation_[guess_index];
    }

    std::array<int, 3> tmp{};
    stripack_trfind_f90(
        idx,
        coords.data(),
        static_cast<int>(num_nodes_),
        xs_.data(), ys_.data(), zs_.data(),
        list_.data(), lptr_.data(), lend_.data(),
        barycentricCoords.data(),
        tmp.data());

    // If all indices are -1 (Fortran index 0), STRIPACK identified all points as coplanar. We can't
    // recover from this, so we ASSERT for this case.
    // If one (the last) index is -1 (Fortran index 0), STRIPACK identified the target point as being
    // outside the convex hull of the source points. We can recover from this by passing the invalid
    // index to the interpolator, which will in turn fill the target result with a missing value.
    constexpr int invalid = -1;
    const size_t numInvalid = std::count(tmp.cbegin(), tmp.cend(), invalid);
    ASSERT(numInvalid == 0 || numInvalid == 1);  // point is enclosed or is outside hull
    if (numInvalid == 1) {
        return false;
    }

    if (reshuffled_) {
      // Invert the randomization, so the returned index matches the coordinates used in construction
      for (size_t i = 0; i < 3; ++i) {
          vertexIndices[i] = inverse_random_permutation_[tmp[i]];
      }
    }
    return true;
}

}  // namespace stripack
