/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Mar 2016

#pragma once

#include <cstddef>

namespace atlas {
namespace grid {
namespace detail {
namespace pl {
namespace classic_gaussian {

/**
 * @brief Compute points per latitude between North pole and equator
 * @param N    [in]  Number of latitudes between pole and equator (Gaussian N
 * number)
 * @param nlon [out] points per latitude
 */
void points_per_latitude_npole_equator(const size_t N, long nlon[]);

/**
 * @brief Compute points per latitude between North pole and equator
 * @param N    [in]  Number of latitudes between pole and equator (Gaussian N
 * number)
 * @param nlon [out] points per latitude
 */
void points_per_latitude_npole_equator(const size_t N, int nlon[]);

/**
 * @brief Compute points per latitude between North pole and South pole
 * @param N    [in]  Number of latitudes between pole and equator (Gaussian N
 * number)
 * @param nlon [out] points per latitude  (size 2*N)
 */
void points_per_latitude_npole_spole(const size_t N, long nlon[]);

/**
 * @brief Compute points per latitude between North pole and South pole
 * @param N    [in]  Number of latitudes between pole and equator (Gaussian N
 * number)
 * @param nlon [out] points per latitude  (size 2*N)
 */
void points_per_latitude_npole_spole(const size_t N, int nlon[]);

}  // namespace classic_gaussian
}  // namespace pl
}  // namespace detail
}  // namespace grid
}  // namespace atlas
