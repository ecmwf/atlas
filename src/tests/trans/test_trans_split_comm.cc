/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */


#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/array.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/output/Gmsh.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/trans/Trans.h"

#include "atlas/util/function/VortexRollup.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::grid::detail::partitioner::TransPartitioner;


namespace atlas {
namespace test {

int getColour(const int world_rank = mpi::comm("world").rank(),
              const int world_size = mpi::comm("world").size()) {
    // 1:3 ratio.
    constexpr float ratio = 0.25;

    if (static_cast<float>(world_rank + 1) / static_cast<float>(world_size) <= ratio) {
        return 0;
    } else {
        return 1;
    }
}

struct AtlasSplitCommEnv : public AtlasTestEnvironment {
    AtlasSplitCommEnv(int argc, char* argv[]): AtlasTestEnvironment(argc, argv) {
        ATLAS_ASSERT(mpi::comm().size() > 1);

        // Split world communicator and set default comm to the split communicator.
        mpi::comm().split(getColour(), "split_comm");
        eckit::mpi::setCommDefault("split_comm");
    }
};


// Each MPI colour in the test does a different resolution transform.
size_t gaussResolutionForComm() {
    static constexpr std::array<size_t, 2> resols = {20, 60};
    return resols[getColour()];
}


CASE("test_trans_split_comm") {
    const size_t resol = gaussResolutionForComm();
    const std::string grid_uid = "O" + std::to_string(resol);

    StructuredGrid g(grid_uid);

    auto N = GaussianGrid(g).N();
    functionspace::Spectral specFS(2 * N - 1);

    functionspace::StructuredColumns gridFS(g,
        grid::Partitioner(new grid::detail::partitioner::TransPartitioner()));

    trans::Trans transIFS(gridFS, specFS);

    Field specField = specFS.createField<double>(option::name("specfield"));

    // Make vortex rollup field on the grid.
    atlas::Field vortexField = gridFS.createField<double>(atlas::option::name("vortex"));
    atlas::Field fieldOut = gridFS.createField<double>(atlas::option::name("out"));
    const auto lonlat = atlas::array::make_view<const double, 2>(vortexField.functionspace().lonlat());
    auto view = atlas::array::make_view<double, 1>(vortexField);

    for (atlas::idx_t ij = 0; ij < view.shape(0); ++ij) {
        constexpr double zz = 1.0;
        view(ij) = atlas::util::function::vortex_rollup(lonlat(ij, 0), lonlat(ij, 1), zz);
    }

    // Go to and from spectral space.
    transIFS.dirtrans(vortexField, specField);
    transIFS.invtrans(specField, fieldOut);

    // Output. Writes two files - one for each split communicator.
    {
        Mesh mesh = StructuredMeshGenerator().generate(g);
        output::Gmsh gmsh(grid_uid + "-grid-split_comm_" + std::to_string(getColour()) + ".msh");
        gmsh.write(mesh);
        gmsh.write(fieldOut, gridFS);
    }
}


}  // namespace test
}  // namespace atlas

//--

int main(int argc, char** argv) {
    return atlas::test::run<atlas::test::AtlasSplitCommEnv>(argc, argv);
}
