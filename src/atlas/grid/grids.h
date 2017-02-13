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


#ifndef atlas_grid_grids_h
#define atlas_grid_grids_h

#include "atlas/grid/Grid.h"
#include "atlas/grid/CustomStructured.h"
#include "atlas/grid/Structured.h"
#include "atlas/grid/Unstructured.h"
#include "atlas/grid/regular/RegularLonLat.h"
#include "atlas/grid/regular/ShiftedLon.h"
#include "atlas/grid/regular/ShiftedLat.h"
#include "atlas/grid/regular/ShiftedLonLat.h"
#include "atlas/grid/regular/RegularGaussian.h"
#include "atlas/grid/reduced/ReducedGaussian.h"
#include "atlas/grid/reduced/ReducedLonLat.h"
#include "atlas/grid/reduced/ClassicGaussian.h"
#include "atlas/grid/reduced/OctahedralGaussian.h"
#include "atlas/grid/reduced/ARPEGE.h"


namespace atlas {
namespace grid {


void load();


void unload();


Grid* grid_from_uid(const std::string& uid);


extern "C" {


    void atlas__grids__load();


}


}  // namespace grid
}  // namespace atlas


#endif
