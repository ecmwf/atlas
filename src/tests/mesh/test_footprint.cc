/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "eckit/log/Bytes.h"

#include "atlas/array/Array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Metadata.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::util;

namespace atlas {
namespace test {

static std::string griduid() {
    return "O64";
}

//-----------------------------------------------------------------------------

CASE("test_footprint") {
    array::ArrayT<double> array(10, 2);
    Log::info() << "array.footprint = " << eckit::Bytes(array.footprint()) << std::endl;

    Field field("field", array::make_datatype<double>(), array::make_shape(10, 2));
    Log::info() << "field.footprint = " << eckit::Bytes(field.footprint()) << std::endl;

    Grid grid(griduid());
    MeshGenerator meshgen("structured");
    Mesh mesh = meshgen.generate(grid);

    Log::info() << "Footprint for mesh generated from grid " << grid.name() << std::endl;
    Log::info() << "mesh.footprint = " << eckit::Bytes(mesh.footprint()) << std::endl;
    Log::info() << "    .nodes.footprint = " << eckit::Bytes(mesh.nodes().footprint()) << std::endl;
    for (idx_t f = 0; f < mesh.nodes().nb_fields(); ++f) {
        Log::info() << "          ." + mesh.nodes().field(f).name() + ".footprint = "
                    << eckit::Bytes(mesh.nodes().field(f).footprint()) << std::endl;
    }

    Log::info() << "    .cells.footprint = " << eckit::Bytes(mesh.cells().footprint()) << std::endl;

    for (idx_t f = 0; f < mesh.cells().nb_fields(); ++f) {
        Log::info() << "          ." + mesh.cells().field(f).name() + ".footprint = "
                    << eckit::Bytes(mesh.cells().field(f).footprint()) << std::endl;
    }

    Log::info() << "          .node_connectivity.footprint = "
                << eckit::Bytes(mesh.cells().node_connectivity().footprint()) << std::endl;
    Log::info() << "          .edge_connectivity.footprint = "
                << eckit::Bytes(mesh.cells().edge_connectivity().footprint()) << std::endl;
    Log::info() << "          .cell_connectivity.footprint = "
                << eckit::Bytes(mesh.cells().cell_connectivity().footprint()) << std::endl;
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
