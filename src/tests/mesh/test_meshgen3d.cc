/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/meshgenerator.h"
#include "atlas/option/Options.h"
#include "atlas/output/Gmsh.h"
#include "tests/AtlasTestEnvironment.h"

using namespace atlas::output;
using namespace atlas::meshgenerator;
using namespace atlas::grid;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

CASE("test_create_mesh_simple") {
    auto mesh = Mesh{Grid{"O32"}};
    Gmsh{"O32.msh"}.write(mesh);
}


CASE("test_create_mesh") {
    Mesh m;

    util::Config opts;
    opts.set("3d", true);            ///< creates links along date-line
    opts.set("include_pole", true);  ///< triangulate the pole point
    StructuredMeshGenerator generate(opts);

    // opts.set("nb_parts",1); // default = 1
    // opts.set("part",    0); // default = 0

    m = generate(Grid("N24"));

    Grid grid = m.grid();
    std::cout << grid.spec() << std::endl;

    Gmsh gmsh("out.msh", util::Config("coordinates", "xyz"));
    gmsh.write(m);

    Field f1("f1", array::make_datatype<double>(), array::make_shape(m.nodes().size()));
    auto v1 = array::make_view<double, 1>(f1);
    for (idx_t i = 0; i < v1.size(); ++i) {
        v1[i] = i + 1;
    }
    gmsh.write(f1);

    std::vector<double> v2(m.nodes().size());
    for (size_t i = 0; i < v2.size(); ++i) {
        v2[i] = i + 1;
    }
    Field f2("f2", v2.data(), array::make_shape(v2.size()));
    gmsh.write(f2);
}

//-----------------------------------------------------------------------------

CASE("METV-3396 / ATLAS-375") {
    auto opts = util::Config("triangulate", true);
    StructuredMeshGenerator generate(opts);
    SECTION("N16 [-180,180]") {
        auto mesh = generate(Grid("N16", RectangularDomain({-180., 180.}, {-90., 90.})));
        Gmsh{"N16.msh"}.write(mesh);
    }
    SECTION("O8 [-45,315]") {
        auto mesh = generate(Grid("O8", RectangularDomain({-45., 315.}, {-90., 90.})));
        Gmsh{"O8.msh"}.write(mesh);
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
