/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>

#include "atlas/array.h"
#include "atlas/grid.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/Locate.h"
#include "atlas/parallel/Collect.h"


#include "tests/AtlasTestEnvironment.h"

//-----------------------------------------------------------------------------


namespace atlas::test {

CASE("real test") {
    const idx_t  remote_index_base = ATLAS_HAVE_FORTRAN;
    const gidx_t global_index_base = 1; // convention
    const int partition_base = 0;

    Log::info() << "remote_index_base = " << remote_index_base << std::endl;

    int mpi_rank = mpi::comm().rank();
    int mpi_size = mpi::comm().size();

    std::vector<gidx_t> receive_gidx;
    std::vector<idx_t>  expected_receive_ridx;
    std::vector<int>    expected_receive_part;

    if (mpi_size == 3) {
        if (mpi_rank == 0) {
            receive_gidx          = {4, 9, 20, 32, 109, 205, 27, 11, 893};
            expected_receive_part = {0, 0, 0,  0,    0,   0,  0,  0,   0};
            expected_receive_ridx = {3, 8, 19, 31, 108, 204, 26, 10, 892};
        }
        if (mpi_rank == 1) {
            receive_gidx          = {7, 52, 305, 1004, 5008};
            expected_receive_part = {0,  0,   0,    0,    2};
            expected_receive_ridx = {6, 51, 304, 1003, 1508};
        }
        if (mpi_rank == 2) {
            receive_gidx          = {989, 2781, 306, 3819, 29};
            expected_receive_part = {0,      1,   0,    2,  0};
            expected_receive_ridx = {988, 1030, 305,  319, 28};
        }
        for( auto& ridx: expected_receive_ridx ) {
            ridx += remote_index_base;
        }
    }
    else {
        Log::warning() << "WARNING: This test needs to be called with mpi_size = 3" << std::endl;
    }

    atlas::Grid grid("O32");
    functionspace::StructuredColumns fs(grid);

    std::vector<idx_t> receive_ridx(receive_gidx.size());
    std::vector<int>   receive_part(receive_gidx.size());

    functionspace::Locator locator(fs);
    locator.locate(
        receive_gidx, global_index_base,
        receive_part, partition_base,
        receive_ridx, remote_index_base);

    EXPECT_EQ(receive_part,expected_receive_part);
    EXPECT_EQ(receive_ridx,expected_receive_ridx);


    // From here, test if we can use this in Collect as a second validation

    parallel::Collect collect;
    collect.setup(receive_part.size(), receive_part.data(),receive_ridx.data(),remote_index_base);

    array::ArrayT<gidx_t> receive(receive_gidx.size());

    collect.execute<gidx_t,1>(fs.global_index(), receive);

    for( int rank=0; rank < mpi_size; ++rank ) {
        mpi::comm().barrier();
        if (rank == mpi_rank) {
            std::cerr << "[" << mpi_rank << "] receive = ";
            receive.dump(std::cerr);
            std::cerr << std::endl;
        }
    }

    {
        auto receive_view = array::make_view<gidx_t,1>(receive);
        for (int i=0; i<receive_gidx.size(); ++i) {
            EXPECT_EQ( receive_view[i], receive_gidx[i]);
        }
    }

}

}  // namespace atlas::test


int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
