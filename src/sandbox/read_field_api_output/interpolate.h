/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/runtime/Exception.h"

#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/interpolation.h"

namespace atlas {

void interpolate(const std::string interpolation_method, const Field& IFS_blocked_f, Field& rad_blocked_f) {
    auto IFS_blocked_fs = IFS_blocked_f.functionspace();
    auto rad_blocked_fs = rad_blocked_f.functionspace();
    auto IFS_grid = IFS_blocked_fs.grid();
    auto rad_grid = rad_blocked_fs.grid();
    if (interpolation_method == "bilinear") {
        auto IFS_nonblocked_fs = functionspace::StructuredColumns{IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs), option::halo(1)};
        auto rad_nonblocked_fs = functionspace::PointCloud{rad_grid, grid::MatchingPartitioner(IFS_blocked_fs) };
        auto interpolation = Interpolation(option::type("structured-bilinear"),  IFS_nonblocked_fs, rad_nonblocked_fs);
        auto IFS_nonblocked_f = IFS_nonblocked_fs.createField(IFS_blocked_f);
        auto rad_nonblocked_f = rad_nonblocked_fs.createField(IFS_blocked_f);
        copy_blocked_to_nonblocked(IFS_blocked_f, IFS_nonblocked_f);
        IFS_nonblocked_f.haloExchange();
        interpolation.execute(IFS_nonblocked_f, rad_nonblocked_f);
        rad_nonblocked_f.haloExchange();
        copy_nonblocked_to_blocked(rad_nonblocked_f, rad_blocked_f);
    }
    else if (interpolation_method == "bicubic") {
        auto IFS_nonblocked_fs = functionspace::StructuredColumns{IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs), option::halo(2)};
        auto rad_nonblocked_fs = functionspace::PointCloud{rad_grid, grid::MatchingPartitioner(IFS_blocked_fs) };
        auto interpolation = Interpolation(option::type("structured-bicubic"),  IFS_nonblocked_fs, rad_nonblocked_fs);
        auto IFS_nonblocked_f = IFS_nonblocked_fs.createField(IFS_blocked_f);
        auto rad_nonblocked_f = rad_nonblocked_fs.createField(IFS_blocked_f);
        copy_blocked_to_nonblocked(IFS_blocked_f, IFS_nonblocked_f);
        IFS_nonblocked_f.haloExchange();
        interpolation.execute(IFS_nonblocked_f, rad_nonblocked_f);
        rad_nonblocked_f.haloExchange();
        copy_nonblocked_to_blocked(rad_nonblocked_f, rad_blocked_f);
    }
    else if (interpolation_method == "cons") {
        auto IFS_mesh = Mesh(IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs));
        auto rad_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_mesh));
        auto IFS_nonblocked_fs = functionspace::NodeColumns{IFS_mesh, option::halo(10)};
        auto rad_nonblocked_fs = functionspace::NodeColumns{rad_mesh, option::halo(1)};
        auto interpolation = Interpolation(option::type("conservative-spherical-polygon"),  IFS_nonblocked_fs, rad_nonblocked_fs);
        auto IFS_nonblocked_f = IFS_nonblocked_fs.createField(IFS_blocked_f);
        auto rad_nonblocked_f = rad_nonblocked_fs.createField(IFS_blocked_f);
        copy_blocked_to_nonblocked(IFS_blocked_f, IFS_nonblocked_f);
        IFS_nonblocked_f.haloExchange();
        interpolation.execute(IFS_nonblocked_f, rad_nonblocked_f);
        rad_nonblocked_f.haloExchange();
        copy_nonblocked_to_blocked(rad_nonblocked_f, rad_blocked_f);
    }
    else if (interpolation_method == "cons2") {
        auto IFS_mesh = Mesh(IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs));
        auto rad_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_mesh));
        auto IFS_nonblocked_fs = functionspace::NodeColumns{IFS_mesh, option::halo(10)};
        auto rad_nonblocked_fs = functionspace::NodeColumns{rad_mesh, option::halo(1)};
        auto interpolation = Interpolation(option::type("conservative-spherical-polygon") | util::Config("order",2),  IFS_nonblocked_fs, rad_nonblocked_fs);
        auto IFS_nonblocked_f = IFS_nonblocked_fs.createField(IFS_blocked_f);
        auto rad_nonblocked_f = rad_nonblocked_fs.createField(IFS_blocked_f);
        copy_blocked_to_nonblocked(IFS_blocked_f, IFS_nonblocked_f);
        IFS_nonblocked_f.haloExchange();
        interpolation.execute(IFS_nonblocked_f, rad_nonblocked_f);
        rad_nonblocked_f.haloExchange();
        copy_nonblocked_to_blocked(rad_nonblocked_f, rad_blocked_f);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

//-----------------------------------------------------------------------------

} // namespace
