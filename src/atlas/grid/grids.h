/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Nov 2014

#ifndef atlas_grids_grids_h
#define atlas_grids_grids_h

#include "atlas/grid/Grid.h"
#include "atlas/grid/global/lonlat/ReducedLonLat.h"
#include "atlas/grid/global/lonlat/RegularLonLat.h"
#include "atlas/grid/global/lonlat/ShiftedLonLat.h"
#include "atlas/grid/global/lonlat/ShiftedLon.h"
#include "atlas/grid/global/lonlat/ShiftedLat.h"
#include "atlas/grid/global/gaussian/Gaussian.h"
#include "atlas/grid/global/gaussian/ReducedGaussian.h"
#include "atlas/grid/global/gaussian/RegularGaussian.h"
#include "atlas/grid/global/gaussian/ClassicGaussian.h"
#include "atlas/grid/global/gaussian/OctahedralGaussian.h"
#include "atlas/grid/global/Structured.h"
#include "atlas/grid/global/CustomStructured.h"
#include "atlas/grid/global/lonlat/LonLat.h"
#include "atlas/grid/global/Unstructured.h"

namespace atlas {
namespace grid {

void load();
void unload();

Grid* grid_from_uid(const std::string& uid);

ReducedGrid* new_reduced_grid(const std::string& identifier);
ReducedGrid* new_reduced_gaussian_grid(const std::vector<int>& nlon);
ReducedGrid* new_gaussian_grid(int N);
ReducedGrid* new_lonlat_grid(int nlon, int nlat);

extern "C"
{
  void atlas__grids__load();
}

} // namespace grid
} // namespace atlas

#endif // atlas_grids_grids_h
