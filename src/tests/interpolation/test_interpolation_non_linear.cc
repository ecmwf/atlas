/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include <algorithm>
#include <limits>

#include "atlas/array.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/NonLinear.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/runtime/Exception.h"

#include "tests/AtlasTestEnvironment.h"


namespace atlas {
namespace test {


const double missingValue    = 42.;
const double missingValueEps = 1e-9;
const double nan             = std::numeric_limits<double>::quiet_NaN();

using field::MissingValue;
using util::Config;

CASE("Interpolation with MissingValue") {
    /*
       Set input field full of 1's, with 9 nodes
         1 ... 1 ... 1
         :     :     :
         1-----m ... 1  m: missing value
         |i   i|     :  i: interpolation on two points, this quadrilateral only
         1-----1 ... 1
     */
    RectangularDomain domain({0, 2}, {0, 2}, "degrees");
    Grid gridA("L90", domain);

    const idx_t nbNodes = 9;
    ATLAS_ASSERT(gridA.size() == nbNodes);

    Mesh meshA = MeshGenerator("structured").generate(gridA);

    functionspace::NodeColumns fsA(meshA);
    Field fieldA = fsA.createField<double>(option::name("A"));

    fieldA.metadata().set("missing_value", missingValue);
    fieldA.metadata().set("missing_value_epsilon", missingValueEps);

    auto viewA = array::make_view<double, 1>(fieldA);
    for (idx_t j = 0; j < fsA.nodes().size(); ++j) {
        viewA(j) = 1;
    }


    // Set output field (2 points)
    functionspace::PointCloud fsB({PointLonLat{0.1, 0.1}, PointLonLat{0.9, 0.9}});


    SECTION("missing-if-all-missing") {
        Interpolation interpolation(Config("type", "finite-element").set("non_linear", "missing-if-all-missing"), fsA,
                                    fsB);

        for (std::string type : {"equals", "approximately-equals", "nan"}) {
            Field fieldB("B", array::make_datatype<double>(), array::make_shape(fsB.size()));
            auto viewB = array::make_view<double, 1>(fieldB);

            fieldA.metadata().set("missing_value_type", type);
            viewA(4) = type == "nan" ? nan : missingValue;

            EXPECT(MissingValue(fieldA));
            interpolation.execute(fieldA, fieldB);

            MissingValue mv(fieldB);
            EXPECT(mv);
            EXPECT(mv(viewB(0)) == false);
            EXPECT(mv(viewB(1)) == false);
        }
    }


    SECTION("missing-if-any-missing") {
        Interpolation interpolation(Config("type", "finite-element").set("non_linear", "missing-if-any-missing"), fsA,
                                    fsB);

        for (std::string type : {"equals", "approximately-equals", "nan"}) {
            Field fieldB("B", array::make_datatype<double>(), array::make_shape(fsB.size()));
            auto viewB = array::make_view<double, 1>(fieldB);

            fieldA.metadata().set("missing_value_type", type);
            viewA(4) = type == "nan" ? nan : missingValue;

            EXPECT(MissingValue(fieldA));
            interpolation.execute(fieldA, fieldB);

            MissingValue mv(fieldB);
            EXPECT(mv);
            EXPECT(mv(viewB(0)));
            EXPECT(mv(viewB(1)));
        }
    }


    SECTION("missing-if-heaviest-missing") {
        Interpolation interpolation(Config("type", "finite-element").set("non_linear", "missing-if-heaviest-missing"),
                                    fsA, fsB);

        for (std::string type : {"equals", "approximately-equals", "nan"}) {
            Field fieldB("B", array::make_datatype<double>(), array::make_shape(fsB.size()));
            auto viewB = array::make_view<double, 1>(fieldB);

            fieldA.metadata().set("missing_value_type", type);
            viewA(4) = type == "nan" ? nan : missingValue;

            EXPECT(MissingValue(fieldA));
            interpolation.execute(fieldA, fieldB);

            MissingValue mv(fieldB);
            EXPECT(mv);
            EXPECT(mv(viewB(0)) == false);
            EXPECT(mv(viewB(1)));
        }
    }
}

CASE("Interpolation of rank 2 field with MissingValue") {

    RectangularDomain domain({0, 2}, {0, 2}, "degrees");
    Grid gridA("L90", domain);

    Mesh meshA = MeshGenerator("structured").generate(gridA);

    functionspace::NodeColumns fsA(meshA);
    Field fieldA = fsA.createField<double>(option::name("A") | option::levels(2) );

    fieldA.metadata().set("missing_value", missingValue);
    fieldA.metadata().set("missing_value_epsilon", missingValueEps);

    auto viewA = array::make_view<double, 2>(fieldA);
    for (idx_t j = 0; j < viewA.shape(0); ++j) {
        for (idx_t k = 0; k < viewA.shape(1); ++k) {
            viewA(j,k) = 1;
        }
    }

    // Set output field (2 points)
    functionspace::PointCloud fsB({PointLonLat{0.1, 0.1}, PointLonLat{0.9, 0.9}});

    SECTION("missing-if-all-missing") {
        Interpolation interpolation(Config("type", "finite-element").set("non_linear", "missing-if-all-missing"), fsA,
                                    fsB);

        for (std::string type : {"equals", "approximately-equals", "nan"}) {
            Field fieldB("B", array::make_datatype<double>(), array::make_shape(fsB.size()));
            auto viewB = array::make_view<double, 1>(fieldB);

            fieldA.metadata().set("missing_value_type", type);
            viewA(4) = type == "nan" ? nan : missingValue;

            EXPECT(MissingValue(fieldA));
            interpolation.execute(fieldA, fieldB);

            MissingValue mv(fieldB);
            EXPECT(mv);
            EXPECT(mv(viewB(0)) == false);
            EXPECT(mv(viewB(1)) == false);

            Field fieldBr2("B", array::make_datatype<double>(), array::make_shape(fsB.size(),2));
            auto viewBr2 = array::make_view<double, 2>(fieldBr2);

            interpolation.execute(fieldA, fieldBr2);

            EXPECT(mv(viewBr2(0,0)) == false);
            EXPECT(mv(viewBr2(1,0)) == false);
            EXPECT(mv(viewBr2(0,1)) == false);
            EXPECT(mv(viewBr2(1,1)) == false);
        }
    }

}


}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
