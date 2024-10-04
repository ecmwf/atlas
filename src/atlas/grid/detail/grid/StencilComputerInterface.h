/*
 * (C) Copyright 2023 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/StencilComputer.h"
#include "atlas/grid/StructuredGrid.h"

extern "C" {
    atlas::grid::ComputeNorth* atlas__grid__ComputeNorth__new(const atlas::StructuredGrid::Implementation* grid, int halo);
    void atlas__grid__ComputeNorth__delete(atlas::grid::ComputeNorth* This);
    atlas::idx_t atlas__grid__ComputeNorth__execute_real32(const atlas::grid::ComputeNorth* This, float y);
    atlas::idx_t atlas__grid__ComputeNorth__execute_real64(const atlas::grid::ComputeNorth* This, double y);

    atlas::grid::ComputeWest* atlas__grid__ComputeWest__new(const atlas::StructuredGrid::Implementation* grid, int halo);
    void atlas__grid__ComputeWest__delete(atlas::grid::ComputeWest* This);
    atlas::idx_t atlas__grid__ComputeWest__execute_real32(const atlas::grid::ComputeWest* This, float x, atlas::idx_t j);
    atlas::idx_t atlas__grid__ComputeWest__execute_real64(const atlas::grid::ComputeWest* This, double x, atlas::idx_t j);

}
