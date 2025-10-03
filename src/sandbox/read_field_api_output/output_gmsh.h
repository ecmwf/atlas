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
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/output/Gmsh.h"
#include "atlas/functionspace/BlockStructuredColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/PointCloud.h"

namespace atlas {

inline void output_gmsh(Field& field, const std::string& gmsh_file, const util::Config& gmsh_config) {
    if (auto blocked_fs = functionspace::BlockStructuredColumns(field.functionspace())) {
        auto mesh = Mesh(blocked_fs.grid(), grid::MatchingPartitioner(blocked_fs));
        auto nonblocked_fs = functionspace::NodeColumns{mesh};
        auto nonblocked_field = nonblocked_fs.createField(field);
        copy_blocked_to_nonblocked(field, nonblocked_field);
        nonblocked_field.haloExchange();
        auto gmsh = output::Gmsh(gmsh_file, gmsh_config);
        gmsh.write(mesh);
        gmsh.write(nonblocked_field);
    }
};

//-----------------------------------------------------------------------------

} // namespace
