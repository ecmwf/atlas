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
#include "atlas/parallel/Collect.h"

#include "tests/AtlasTestEnvironment.h"


/// POD: Type to test
using POD = double;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

template <typename T, size_t N>
std::vector<T> vec(const T (&list)[N]) {
    return std::vector<T>(list, list + N);
}


struct Fixture {
    /*
       
    Data:
    [0] : 9[2,0], 1[0,1], 2[0,2], 3[0,3], 4[1,1], 20[1,2]
    [1] : 3[0,3], 4[1,1], 5[1,2], 6[1,3], 7[2,1],  8[2,3]
    [2] : 5[1,2], 6[1,3], 7[2,2], 8[2,3], 9[2,4],  1[0,1], 2[0,2]


                        1[0,1]     2[0,2]

                9[2,4]                     3[0,3]

                8[2,3]                    4[1,1]

                        7[2,2]          5[1,2]

                            6[1,3]


    Receive gidx:
    [0] : 7[2,2], 5[1,2], 6[1,3]
    [1] : 9[2,4], 3[0,3], 8[2,3], 4[1,1]
    [2] : 1[0,1], 2[0,2]
 
    */
    Fixture() {
        rank           = static_cast<int>(mpi::comm().rank());
        comm_size      = static_cast<int>(mpi::comm().size());
        int nnodes_c[] = {6, 6, 7};
        nb_nodes       = vec(nnodes_c);
        Nl             = nb_nodes[rank];
        int ridx_base  = 0;
        switch (mpi::comm().rank()) {
            case 0: {  //./----> extra ghost point with nonstandard gidx
                int part_c[]    = {2, 0, 0, 0, 1, 2};
                part            = vec(part_c);
                idx_t ridx_c[]  = {4, 1, 2, 3, 1, 3};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {9, 1, 2, 3, 4, 20};
                gidx            = vec(gidx_c);
                int recv_part_c[] = {2, 1, 1};
                recv_part       = vec(recv_part_c);
                idx_t recv_ridx_c[] = {2, 2, 3};
                recv_ridx       = vec(recv_ridx_c);
                break;
            }
            case 1: {
                int part_c[]    = {0, 1, 1, 1, 2, 2};
                part            = vec(part_c);
                idx_t ridx_c[]  = {3, 1, 2, 3, 2, 3};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {3, 4, 5, 6, 7, 8};
                gidx            = vec(gidx_c);
                int recv_part_c[] = {2, 0, 2, 1};
                recv_part       = vec(recv_part_c);
                idx_t recv_ridx_c[] = {4, 3, 3, 1};
                recv_ridx       = vec(recv_ridx_c);
                break;
            }
            case 2: {
                int part_c[]    = {1, 1, 2, 2, 2, 0, 0};
                part            = vec(part_c);
                idx_t ridx_c[]  = {2, 3, 2, 3, 4, 1, 2};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {5, 6, 7, 8, 9, 1, 2};
                gidx            = vec(gidx_c);
                int recv_part_c[] = {0, 0};
                recv_part       = vec(recv_part_c);
                idx_t recv_ridx_c[] = {1, 2};
                recv_ridx       = vec(recv_ridx_c);
                break;
            }
        }
        recv_size = recv_part.size();
        collect.setup(recv_size, recv_part.data(), recv_ridx.data(), ridx_base);
    }
    parallel::Collect collect;
    std::vector<int> nb_nodes;
    std::vector<int> part;
    std::vector<idx_t> ridx;
    std::vector<gidx_t> gidx;

    std::vector<int>    recv_part;
    std::vector<idx_t>  recv_ridx;

    int recv_size;
    int ridx_base;
    int Nl;
    int root;
    int rank;
    int comm_size;
};

//-----------------------------------------------------------------------------

CASE("test_all_to_all rank 1") {
    Fixture f;

    array::ArrayT<POD> arr_loc(f.Nl);
    array::ArrayT<POD> arr_recv(f.recv_size);

    {
        auto loc = array::make_view<POD,1>(arr_loc);
        for (int j = 0; j < f.Nl; ++j) {
            loc[j] = f.gidx[j];
        }
    }

    f.collect.execute<POD,1>(arr_loc, arr_recv);

    {
        auto recv = array::make_view<POD,1>(arr_recv);
        if (f.rank == 0) {
            EXPECT_EQ(recv(0), 7);
            EXPECT_EQ(recv(1), 5);
            EXPECT_EQ(recv(2), 6);
        }
        if (f.rank == 1) {
            EXPECT_EQ(recv(0), 9);
            EXPECT_EQ(recv(1), 3);
            EXPECT_EQ(recv(2), 8);
            EXPECT_EQ(recv(3), 4);
        }
        if (f.rank == 2) {
            EXPECT_EQ(recv(0), 1);
            EXPECT_EQ(recv(1), 2);
        }
    }
}

CASE("test_all_to_all rank 2") {
    Fixture f;

    array::ArrayT<POD> arr_loc(f.Nl, 2);
    array::ArrayT<POD> arr_recv(f.recv_size, 2);

    {
        auto loc = array::make_view<POD,2>(arr_loc);
        for (int j = 0; j < f.Nl; ++j) {
            loc(j,0) = f.gidx[j];
            loc(j,1) = 10 * f.gidx[j];
        }
    }

    f.collect.execute<POD,2>(arr_loc, arr_recv);

    {
        auto recv = array::make_view<POD,2>(arr_recv);
        if (f.rank == 0) {
            EXPECT_EQ(recv(0,0), 7); EXPECT_EQ(recv(0,1), 70);
            EXPECT_EQ(recv(1,0), 5); EXPECT_EQ(recv(1,1), 50);
            EXPECT_EQ(recv(2,0), 6); EXPECT_EQ(recv(2,1), 60);
        }
        if (f.rank == 1) {
            EXPECT_EQ(recv(0,0), 9); EXPECT_EQ(recv(0,1), 90);
            EXPECT_EQ(recv(1,0), 3); EXPECT_EQ(recv(1,1), 30);
            EXPECT_EQ(recv(2,0), 8); EXPECT_EQ(recv(2,1), 80);
            EXPECT_EQ(recv(3,0), 4); EXPECT_EQ(recv(3,1), 40);
        }
        if (f.rank == 2) {
            EXPECT_EQ(recv(0,0), 1); EXPECT_EQ(recv(0,1), 10);
            EXPECT_EQ(recv(1,0), 2); EXPECT_EQ(recv(1,1), 20);
        }

    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/util/Locate.h"

namespace atlas::test {

CASE("real test") {
    int mpi_rank = mpi::comm().rank();
    int mpi_size = mpi::comm().size();

    std::vector<gidx_t> receive_gidx;
    if (mpi_rank == 0) {
        receive_gidx = {4, 9, 20, 32, 109, 205, 27, 11, 893};
    }
    if (mpi_rank == 1) {
        receive_gidx = {7, 52, 305, 1004, 5008};
    }
    if (mpi_rank == 2) {
        receive_gidx = {989, 2781, 306, 3819, 29};
    }

    atlas::Grid grid("O32");
    functionspace::StructuredColumns fs(grid);

    std::vector<idx_t> receive_ridx;
    std::vector<int> receive_part;

    atlas::util::locate(fs, receive_gidx, receive_part, receive_ridx, 0);


    parallel::Collect collect;
    collect.setup(receive_part.size(), receive_part.data(),receive_ridx.data(),0);

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
