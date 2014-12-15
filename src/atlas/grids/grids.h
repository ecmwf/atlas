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
/// @date Nov 2014

#ifndef atlas_grids_grids_h
#define atlas_grids_grids_h

#include "atlas/grids/GaussianGrid.h"
#include "atlas/grids/LonLatGrid.h"
#include "atlas/grids/ReducedGaussianGrid.h"
#include "atlas/grids/ReducedGrid.h"
#include "atlas/grids/ReducedLonLatGrid.h"
#include "atlas/grids/Unstructured.h"
#include "atlas/grids/rgg/rgg.h"

/// @todo remaining to be reworked
#include "atlas/grids/PolarStereoGraphic.h"
#include "atlas/grids/RotatedLatLon.h"


namespace atlas {
namespace grids {

void load();
void unload();

} // namespace grids
} // namespace atlas

extern "C"
{
  void atlas__grids__load();
}

#endif // atlas_grids_grids_h
