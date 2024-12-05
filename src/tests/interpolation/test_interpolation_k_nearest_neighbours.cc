/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#include <cmath>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/MD5.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Config.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::NodeColumns;
using atlas::util::Config;

namespace atlas {
namespace test {

class Access {
public:
    Access(const Interpolation& interpolation): interpolation_{interpolation} {}
    const Interpolation::Implementation::Matrix& matrix() const { return *interpolation_.get()->matrix_; }
    Interpolation interpolation_;

    std::string hash() {
        eckit::MD5 hash;
        const auto m      = atlas::linalg::make_host_view<eckit::linalg::Scalar,eckit::linalg::Index>(matrix());
        const auto outer  = m.outer();
        const auto index  = m.inner();
        const auto weight = m.value();

        idx_t rows = static_cast<idx_t>(m.rows());
        hash.add(rows);
        for (idx_t r = 0; r < rows; ++r) {
            //v_tgt( r ) = 0.;
            for (idx_t c = outer[r]; c < outer[r + 1]; ++c) {
                hash.add(c);
                idx_t n = index[c];
                hash.add(n);
                //Value w = static_cast<Value>( weight[c] );
                //v_tgt( r ) += w * v_src( n );
            }
        }
        return hash.digest();
    }
};


CASE("test_interpolation_k_nearest_neighbours") {
    Log::info().precision(16);

    Grid gridA("O32");
    Grid gridB("O64");
    Grid gridC("O32");

    Interpolation a(option::type("k-nearest-neighbours"), gridA, gridB);
    Interpolation b(option::type("k-nearest-neighbours"), gridB, gridC);
    ATLAS_DEBUG_VAR(Access{a}.hash());
    ATLAS_DEBUG_VAR(Access{b}.hash());

    // Following EXPECT are commented because they are not bit-identical across platforms.
    // This indicates that not every platform finds the same nearest neighbours!!!
    //EXPECT( Access{a}.hash() == "5ecc9d615dcf7112b5a97b38341099a5" );
    //EXPECT( Access{b}.hash() == "bf4b4214ef50387ba00a1950b95e0d93" );
}

//-----------------------------------------------------------------------------

CASE("test_multiple_fs") {
    Grid grid1("L90x45");
    Grid grid2("O8");

    Mesh mesh1 = StructuredMeshGenerator().generate(grid1);
    Mesh mesh2 = StructuredMeshGenerator().generate(grid2);

    functionspace::NodeColumns fs11(mesh1, option::halo(1));
    functionspace::NodeColumns fs12(mesh1, option::halo(2));

    auto fs1 = fs11;
    functionspace::NodeColumns fs2(mesh2, option::halo(1));

    Interpolation interpolation12(Config("type", "k-nearest-neighbours") | Config("k-nearest-neighbours", 5), fs1, fs2);

    auto f1 = fs1.createField<double>(Config("name", "source"));
    auto f2 = fs2.createField<double>(Config("name", "target"));

    auto v1 = array::make_view<double, 1>(f1);
    v1.assign(1.);

    interpolation12.execute(f1, f2);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
