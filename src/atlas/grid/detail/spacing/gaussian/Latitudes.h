/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @file   Latitudes.h
/// @author Willem Deconinck
/// @date   Jan 2014

#pragma once

#include <cstddef>

namespace atlas {
namespace grid {
namespace spacing {
namespace gaussian {

/**
 * @brief Compute gaussian latitudes between North pole and equator
 * @param N         [in]  Number of latitudes between pole and equator (Gaussian
 * N number)
 * @param latitudes [out] latitudes in degrees
 */
void gaussian_latitudes_npole_equator(const size_t N, double latitudes[]);

/**
 * @brief Compute gaussian latitudes and quadrature weights between North pole
 * and equator
 * @param N         [in]  Number of latitudes between pole and equator (Gaussian
 * N number)
 * @param latitudes [out] latitudes in degrees
 * @param weights   [out] quadrature weights
 */
void gaussian_quadrature_npole_equator(const size_t N, double latitudes[], double weights[]);

/**
 * @brief Compute gaussian latitudes between North pole and South pole
 * @param N         [in]  Number of latitudes between pole and equator (Gaussian
 * N number)
 * @param latitudes [out] latitudes in degrees  (size 2*N)
 */
void gaussian_latitudes_npole_spole(const size_t N, double latitudes[]);

/**
 * @brief Compute gaussian latitudes and quadrature weights between North pole
 * and South pole
 * @param N         [in]  Number of latitudes between pole and equator (Gaussian
 * N number)
 * @param latitudes [out] latitudes in degrees (size 2*N)
 * @param weights   [out] quadrature weights   (size 2*N)
 */
void gaussian_quadrature_npole_spole(const size_t N, double latitudes[], double weights[]);

}  // namespace gaussian
}  // namespace spacing
}  // namespace grid
}  // namespace atlas
