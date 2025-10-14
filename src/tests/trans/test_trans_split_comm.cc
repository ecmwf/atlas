/*
 * (C) Crown Copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */


#include <chrono>
#include <iostream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "eckit/mpi/Comm.h"
#include "eckit/mpi/Parallel.h"

#include "atlas/grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/util/Config.h"

#if ATLAS_HAVE_TRANS
#include "atlas/library/config.h"
#include "atlas/trans/ifs/TransIFS.h"
#include "atlas/trans/ifs/TransIFSNodeColumns.h"
#include "atlas/trans/ifs/TransIFSStructuredColumns.h"
#if ATLAS_HAVE_ECTRANS
#include "ectrans/transi.h"
#else
#include "transi/trans.h"
#endif
#endif

#include "tests/AtlasTestEnvironment.h"


using atlas::grid::detail::partitioner::TransPartitioner;


namespace atlas {
namespace test {

int getColour() {
  return atlas::mpi::comm().rank() % 2;
}

struct AtlasTransEnvironment : public AtlasTestEnvironment {
    AtlasTransEnvironment(int argc, char* argv[]): AtlasTestEnvironment(argc, argv) {
        if (eckit::mpi::comm().size() == 1) {
            trans_use_mpi(false);
            // TODO(JC): Abort
        }

        // Split world communicator and set Atlas default to the split.
        atlas::mpi::comm().split(getColour(), "split_comm");
        eckit::mpi::setCommDefault("split_comm");
    }

    ~AtlasTransEnvironment() { trans_finalize(); }

};


CASE("test_trans_split_comm") {

    const atlas::mpi::Comm& comm = atlas::mpi::comm();
    const std::size_t no_mpi_ranks = comm.size();
    std::size_t mpi_rank = comm.rank();

    std::cout << "communicator: " << comm.name() << std::endl;
    std::cout << "no. of MPI ranks: " << no_mpi_ranks << std::endl;
    std::cout << "MPI rank: " << mpi_rank << std::endl;

    // list of MPI communicators
    std::vector<std::string> ds_comms = eckit::mpi::listComms();
    std::cout << ds_comms << std::endl;

    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    //--

    EXPECT(grid::Partitioner::exists("ectrans"));

    Grid g("N80");

    EXPECT(StructuredGrid(g).ny() == 160);

    auto trans_partitioner = new TransPartitioner();
    grid::Partitioner partitioner(trans_partitioner);
    grid::Distribution distribution(g, partitioner);

    trans::TransIFS trans(g, 159);
    ::Trans_t* t = trans;

    ATLAS_DEBUG_VAR(trans.truncation());
    EXPECT(trans.truncation() == 159);

    EXPECT(t->nproc == int(atlas::mpi::comm().size()));
    EXPECT(t->myproc == int(atlas::mpi::comm().rank() + 1));

    // all tasks do the same, so only one needs to check
    if (atlas::mpi::comm().rank() == 0) {

        int max_nb_regions_EW(0);
        for (int j = 0; j < trans_partitioner->nb_bands(); ++j) {
            max_nb_regions_EW = std::max(max_nb_regions_EW, trans_partitioner->nb_regions(j));
        }

        EXPECT(t->n_regions_NS == trans_partitioner->nb_bands());
        EXPECT(t->n_regions_EW == max_nb_regions_EW);

        EXPECT(distribution.nb_partitions() == idx_t(atlas::mpi::comm().size()));
        EXPECT(idx_t(distribution.size()) == g.size());

        std::vector<int> npts(distribution.nb_partitions(), 0);

        for (idx_t j = 0; j < g.size(); ++j) {
            ++npts[distribution.partition(j)];
        }

        EXPECT(t->ngptotg == g.size());
        EXPECT(t->ngptot == npts[atlas::mpi::comm().rank()]);
        EXPECT(t->ngptotmx == *std::max_element(npts.begin(), npts.end()));

        for (int j = 0; j < trans_partitioner->nb_bands(); ++j) {
            EXPECT(t->n_regions[j] == trans_partitioner->nb_regions(j));
        }
    }

}



}  // namespace test
}  // namespace atlas

//--

int main(int argc, char** argv) {
    return atlas::test::run<atlas::test::AtlasTransEnvironment>(argc, argv);
}
