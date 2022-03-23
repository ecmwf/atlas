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
            if (fieldB.metadata().has("missing_value"))
                std::cout <<  "viewB MissingValue value: " << fieldB.metadata().getDouble("missing_value") << std::endl;
            std::cout <<  "mv(viewB(0)) " << mv(viewB(0));
            std::cout <<  " viewB(0) " << viewB(0) << std::endl;
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

    int nlevels = 3;
    functionspace::NodeColumns fsA(meshA);
    Field fieldA = fsA.createField<double>(option::name("A") | option::levels(nlevels) );

    fieldA.metadata().set("missing_value", missingValue);
    fieldA.metadata().set("missing_value_epsilon", missingValueEps);

    auto viewA = array::make_view<double, 2>(fieldA);
    for (idx_t j = 0; j < viewA.shape(0); ++j) {
        viewA(j,0) = 10;
        viewA(j,1) = missingValue;
        viewA(j,2) = 10;
        //viewA(j,0) = 10+j;
        //viewA(j,1) = 20+j;
        //viewA(j,2) = 30+j;
        std::cout <<  "viewA(" << j << ",0) " << viewA(j,0);
        std::cout << " viewA(" << j << ",1) " << viewA(j,1);
        std::cout << " viewA(" << j << ",2) " << viewA(j,2) << std::endl;
    }

    const array::ArraySpec spec(array::ArrayShape{fieldA.shape(0)}, array::ArrayStrides{fieldA.shape(1)});

    // Set output field (2 points)
    functionspace::PointCloud fsB({PointLonLat{0.1, 0.1}, PointLonLat{0.9, 0.9}});

    SECTION("check wrapped data can be indexed with []") {
        for (int lev=0; lev < nlevels; ++lev) {
            double* data = const_cast<double*>(fieldA.array().data<double>()) + lev;
            atlas::Field fieldA_slice = atlas::Field("s", data, spec);
            fieldA_slice.metadata().set("missing_value", fieldA.metadata().get<double>("missing_value"));

            auto FlatViewA1 = array::make_view<double, 1>(fieldA_slice);
            auto viewA_slice = viewA.slice(array::Range::all(), lev);
            EXPECT(FlatViewA1.shape(0) == viewA_slice.shape(0));
            for (idx_t j = 0; j < FlatViewA1.shape(0); ++j) {
                EXPECT(FlatViewA1[j] == viewA_slice(j));
            }
        }
    }

    SECTION("missing-if-all-missing") {
        Interpolation interpolation(Config("type", "finite-element").set("non_linear", "missing-if-all-missing"), fsA,
                                    fsB);

        for (std::string type : {"equals", "approximately-equals", "nan"}) {

            fieldA.metadata().set("missing_value_type", type);
            for (idx_t j = 0; j < viewA.shape(0); ++j) {
                viewA(j,1) = type == "nan" ? nan : missingValue;
            }

            EXPECT(MissingValue(fieldA));

            Field fieldB = fsB.createField<double>(option::name("B") | option::levels(nlevels) );
            auto viewB = array::make_view<double, 2>(fieldB);

            interpolation.execute(fieldA, fieldB);

            MissingValue mv(fieldB);
            if (fieldB.metadata().has("missing_value"))
                std::cout <<  "viewB MissingValue value: " << fieldB.metadata().getDouble("missing_value") << std::endl;
            for (idx_t j = 0; j < viewB.shape(0); ++j) {
                std::cout <<  "viewB(" << j << ",0) " << viewB(j,0);
                std::cout << " viewB(" << j << ",1) " << viewB(j,1);
                std::cout << " viewB(" << j << ",2) " << viewB(j,2) << std::endl;
            }
            EXPECT(mv(viewB(0,0)) == false);
            EXPECT(mv(viewB(1,0)) == false);
            EXPECT(mv(viewB(0,1)));
            EXPECT(mv(viewB(1,1)));
            EXPECT(mv(viewB(0,2)) == false);
            EXPECT(mv(viewB(1,2)) == false);
        }
    }

}


}  // namespace test
}  // namespace atlas


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
