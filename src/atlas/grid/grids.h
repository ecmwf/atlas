/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date Nov 2014


#ifndef atlas_grid_grids_h
#define atlas_grid_grids_h

#include "atlas/grid/Grid.h"
#include "atlas/grid/CustomStructured.h"
#include "atlas/grid/Structured.h"
#include "atlas/grid/Unstructured.h"
#include "atlas/grid/gaussian/ClassicGaussian.h"
#include "atlas/grid/gaussian/Gaussian.h"
#include "atlas/grid/gaussian/OctahedralGaussian.h"
#include "atlas/grid/gaussian/ReducedGaussian.h"
#include "atlas/grid/gaussian/RegularGaussian.h"
#include "atlas/grid/lonlat/LonLat.h"
#include "atlas/grid/lonlat/ReducedLonLat.h"
#include "atlas/grid/lonlat/RegularLonLat.h"
#include "atlas/grid/lonlat/ShiftedLat.h"
#include "atlas/grid/lonlat/ShiftedLon.h"
#include "atlas/grid/lonlat/ShiftedLonLat.h"


namespace atlas {
namespace grid {

Grid* grid_from_uid(const std::string& uid);

}  // namespace grid
}  // namespace atlas


#endif
