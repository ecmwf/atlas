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

void interpolate(Interpolation interpolation, bool on_device, const FieldSet& IFS_blocked_fields, FieldSet& rad_blocked_fields) {
    FieldSet IFS_nonblocked_fields;
    FieldSet rad_nonblocked_fields;
    auto IFS_nonblocked_fs = interpolation.source();
    auto rad_nonblocked_fs = interpolation.target();
    ATLAS_TRACE_SCOPE("create nonblocked fields "+std::string(on_device?"[device]":"[host]")){
        for(int i=0; i<IFS_blocked_fields.size(); ++i) {
            IFS_nonblocked_fields.add(IFS_nonblocked_fs.createField(IFS_blocked_fields[i]));
            rad_nonblocked_fields.add(rad_nonblocked_fs.createField(IFS_blocked_fields[i]));
        }
        if (on_device) {
            IFS_nonblocked_fields.allocateDevice();
            rad_nonblocked_fields.allocateDevice();
        }
    }
    ATLAS_TRACE_SCOPE("copy_blocked_to_nonblocked "+std::string(on_device?"[device]":"[host]")){
        copy_blocked_to_nonblocked(IFS_blocked_fields, IFS_nonblocked_fields, on_device);
        IFS_nonblocked_fields.haloExchange(on_device);
    }
    ATLAS_TRACE_SCOPE("nonblocked interpolation "+std::string(on_device?"[device]":"[host]")) {
        interpolation.execute(IFS_nonblocked_fields, rad_nonblocked_fields);
    }
    ATLAS_TRACE_SCOPE("copy_nonblocked_to_blocked "+std::string(on_device?"[device]":"[host]")){
        copy_nonblocked_to_blocked(rad_nonblocked_fields, rad_blocked_fields, on_device);
    }
}

Interpolation create_interpolation(const std::string& interpolation_method, bool on_device, FunctionSpace IFS_blocked_fs, FunctionSpace rad_blocked_fs) {
    auto IFS_grid = IFS_blocked_fs.grid();
    auto rad_grid = rad_blocked_fs.grid();
    if (interpolation_method == "bilinear") {
        auto IFS_nonblocked_fs = functionspace::StructuredColumns{IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs), option::halo(1)};
        auto rad_nonblocked_fs = functionspace::PointCloud{rad_grid, grid::MatchingPartitioner(IFS_blocked_fs) };
        return Interpolation(option::type("structured-bilinear"),  IFS_nonblocked_fs, rad_nonblocked_fs);
    }
    else if (interpolation_method == "bicubic") {
        auto IFS_nonblocked_fs = functionspace::StructuredColumns{IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs), option::halo(2)};
        auto rad_nonblocked_fs = functionspace::PointCloud{rad_grid, grid::MatchingPartitioner(IFS_blocked_fs) };
        return Interpolation(option::type("structured-bicubic"),  IFS_nonblocked_fs, rad_nonblocked_fs);
    }
    else if (interpolation_method == "cons") {
        auto IFS_mesh = Mesh(IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs));
        auto rad_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_mesh));
        auto IFS_nonblocked_fs = functionspace::NodeColumns{IFS_mesh, option::halo(10)};
        auto rad_nonblocked_fs = functionspace::NodeColumns{rad_mesh, option::halo(1)};
        return Interpolation(option::type("conservative-spherical-polygon"),  IFS_nonblocked_fs, rad_nonblocked_fs);
    }
    else if (interpolation_method == "cons2") {
        auto IFS_mesh = Mesh(IFS_grid, grid::MatchingPartitioner(IFS_blocked_fs));
        auto rad_mesh = Mesh(rad_grid, grid::MatchingPartitioner(IFS_mesh));
        auto IFS_nonblocked_fs = functionspace::NodeColumns{IFS_mesh, option::halo(10)};
        auto rad_nonblocked_fs = functionspace::NodeColumns{rad_mesh, option::halo(1)};
        return Interpolation(option::type("conservative-spherical-polygon") | util::Config("order",2),  IFS_nonblocked_fs, rad_nonblocked_fs);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

void interpolate(const std::string interpolation_method, bool on_device, const Field& IFS_blocked_f, Field& rad_blocked_f) {
    auto IFS_blocked_fs = IFS_blocked_f.functionspace();
    auto rad_blocked_fs = rad_blocked_f.functionspace();

    auto interpolation = create_interpolation(interpolation_method, on_device, IFS_blocked_fs, rad_blocked_fs);
    {
        FieldSet IFS_nonblocked_fields;
        FieldSet rad_nonblocked_fields;
        auto IFS_nonblocked_fs = interpolation.source();
        auto rad_nonblocked_fs = interpolation.target();
        auto IFS_nonblocked_f = IFS_nonblocked_fs.createField(IFS_blocked_f);
        auto rad_nonblocked_f = rad_nonblocked_fs.createField(IFS_blocked_f);
        copy_blocked_to_nonblocked(IFS_blocked_f, IFS_nonblocked_f, on_device);
        IFS_nonblocked_fields.haloExchange(on_device);
        interpolation.execute(IFS_nonblocked_f, rad_nonblocked_f);
        rad_nonblocked_f.haloExchange(on_device);
        copy_nonblocked_to_blocked(rad_nonblocked_f, rad_blocked_f, on_device);
    };
}

void interpolate(const std::string interpolation_method, bool on_device, const FieldSet& IFS_blocked_fields, FieldSet& rad_blocked_fields) {
    auto IFS_blocked_fs = IFS_blocked_fields[0].functionspace();
    auto rad_blocked_fs = rad_blocked_fields[0].functionspace();

    auto interpolation = create_interpolation(interpolation_method, on_device, IFS_blocked_fs, rad_blocked_fs);
    {
        FieldSet IFS_nonblocked_fields;
        FieldSet rad_nonblocked_fields;
        auto IFS_nonblocked_fs = interpolation.source();
        auto rad_nonblocked_fs = interpolation.target();
        for(int i=0; i<IFS_blocked_fields.size(); ++i) {
            IFS_nonblocked_fields.add(IFS_nonblocked_fs.createField(IFS_blocked_fields[i]));
            rad_nonblocked_fields.add(rad_nonblocked_fs.createField(IFS_blocked_fields[i]));
        }
        copy_blocked_to_nonblocked(IFS_blocked_fields, IFS_nonblocked_fields, on_device);
        IFS_nonblocked_fields.haloExchange(on_device);
        interpolation.execute(IFS_nonblocked_fields, rad_nonblocked_fields);
        rad_nonblocked_fields.haloExchange(on_device);
        copy_nonblocked_to_blocked(rad_nonblocked_fields, rad_blocked_fields, on_device);
    }
}

//-----------------------------------------------------------------------------

} // namespace
