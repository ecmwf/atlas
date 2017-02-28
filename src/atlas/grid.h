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

#pragma once

#include "atlas/grid/Grid.h"

namespace atlas {
namespace grid {


void load();


void unload();


detail::grid::Grid* grid_from_uid(const std::string& uid);


extern "C" {

    void atlas__grids__load();

}


}  // namespace grid
}  // namespace atlas
