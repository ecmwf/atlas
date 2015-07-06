/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Jan 2014

#ifndef atlas_grids_GaussianLatitudes_h
#define atlas_grids_GaussianLatitudes_h

#include <cstddef>

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

void gaussian_latitudes_npole_equator (const size_t N, double[]);
void gaussian_latitudes_npole_spole   (const size_t N, double[]);

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif // atlas_grids_GaussianLatitudes_h
