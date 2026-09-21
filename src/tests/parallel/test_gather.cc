/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <algorithm>
#include <cmath>
#include <sstream>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/library/config.h"
#include "atlas/parallel/GatherScatter.h"
#include "atlas/parallel/mpi/mpi.h"
#include "eckit/system/ResourceUsage.h"
#include "eckit/utils/Translator.h"

#include "tests/AtlasTestEnvironment.h"

/// POD: Type to test
using POD = double;

namespace atlas {
namespace test {

constexpr size_t MB = 1024 * 1024;
constexpr size_t GB = 1024 * MB;

size_t peakMemory() {
    return eckit::system::ResourceUsage().maxResidentSetSize();
}

//-----------------------------------------------------------------------------

template <typename T, size_t N>
std::vector<T> vec(const T (&list)[N]) {
    return std::vector<T>(list, list + N);
}

struct Fixture {
    Fixture() {
        rank           = static_cast<int>(mpi::comm().rank());
        comm_size      = static_cast<int>(mpi::comm().size());
        int nnodes_c[] = {6, 6, 7};
        nb_nodes       = vec(nnodes_c);
        Nl             = nb_nodes[rank];
        switch (mpi::comm().rank()) {
            case 0: {  //./----> extra ghost point with nonstandard gidx
                int part_c[]    = {2, 0, 0, 0, 1, 2};
                part            = vec(part_c);
                idx_t ridx_c[]  = {4, 1, 2, 3, 1, 3};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {9, 1, 2, 3, 4, 20};
                gidx            = vec(gidx_c);
                break;
            }
            case 1: {
                int part_c[]    = {0, 1, 1, 1, 2, 2};
                part            = vec(part_c);
                idx_t ridx_c[]  = {3, 1, 2, 3, 2, 3};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {3, 4, 5, 6, 7, 8};
                gidx            = vec(gidx_c);
                break;
            }
            case 2: {
                int part_c[]    = {1, 1, 2, 2, 2, 0, 0};
                part            = vec(part_c);
                idx_t ridx_c[]  = {2, 3, 2, 3, 4, 1, 2};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {5, 6, 7, 8, 9, 1, 2};
                gidx            = vec(gidx_c);
                break;
            }
        }
        gather_scatter.setup(part.data(), ridx.data(), 0, gidx.data(), Nl);
    }
    parallel::GatherScatter gather_scatter;
    std::vector<int> nb_nodes;
    std::vector<int> part;
    std::vector<idx_t> ridx;
    std::vector<gidx_t> gidx;

    int Nl;
    int root;
    int rank;
    int comm_size;

    int Ng() { return rank == root ? gather_scatter.glb_dof() : 0; }
};

struct SparseGlobalIndexFixture {
    SparseGlobalIndexFixture() {
        rank      = static_cast<int>(mpi::comm().rank());
        comm_size = static_cast<int>(mpi::comm().size());
        int nnodes_c[] = {4, 4, 4};
        nb_nodes       = vec(nnodes_c);
        Nl             = nb_nodes[rank];
        switch (rank) {
            case 0: {
                int part_c[]    = {0, 0, 1, 2};
                part            = vec(part_c);
                idx_t ridx_c[]  = {0, 1, 0, 0};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {10, 30, 20, 40};
                gidx            = vec(gidx_c);
                break;
            }
            case 1: {
                int part_c[]    = {1, 1, 0, 2};
                part            = vec(part_c);
                idx_t ridx_c[]  = {0, 1, 1, 1};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {20, 50, 30, 1439999999};
                gidx            = vec(gidx_c);
                break;
            }
            case 2: {
                int part_c[]    = {2, 2, 0, 1};
                part            = vec(part_c);
                idx_t ridx_c[]  = {0, 1, 0, 1};
                ridx            = vec(ridx_c);
                gidx_t gidx_c[] = {40, 1439999999, 10, 50};
                gidx            = vec(gidx_c);
                break;
            }
        }
        gather_scatter.setup(part.data(), ridx.data(), 0, gidx.data(), Nl);
    }

    parallel::GatherScatter gather_scatter;
    std::vector<int> nb_nodes;
    std::vector<int> part;
    std::vector<idx_t> ridx;
    std::vector<gidx_t> gidx;

    int Nl;
    int root;
    int rank;
    int comm_size;

    int Ng() { return rank == root ? gather_scatter.glb_dof() : 0; }
};

struct SetupOptionsFixture {
    SetupOptionsFixture() {
        rank = static_cast<int>(mpi::comm().rank());
        EXPECT_EQ(mpi::comm().size(), 3);

        switch (rank) {
            case 0:
                part        = {0, 0, 1, 2};
                remote_idx  = {0, 1, 0, 0};
                dense_gidx  = {1, 4, 2, 3};
                break;
            case 1:
                part        = {1, 1, 0, 2};
                remote_idx  = {0, 1, 0, 0};
                dense_gidx  = {2, 5, 1, 3};
                break;
            case 2:
                part        = {2, 2, 0, 1};
                remote_idx  = {0, 1, 0, 0};
                dense_gidx  = {3, 6, 1, 2};
                break;
        }

        mask = {0, 0, 1, 1};
        sparse_gidx.reserve(dense_gidx.size());
        constexpr gidx_t sparse_indices[] = {0, 2, 4, 7, 9, 11, 13};
        for (const gidx_t gidx : dense_gidx) {
            sparse_gidx.emplace_back(sparse_indices[gidx]);
        }
    }

    int rank;
    std::vector<int> part;
    std::vector<idx_t> remote_idx;
    std::vector<gidx_t> dense_gidx;
    std::vector<gidx_t> sparse_gidx;
    std::vector<int> mask;
};

void test_mapping(parallel::GatherScatter& gather_scatter, const SetupOptionsFixture& fixture,
                  const std::vector<gidx_t>& global_index, const std::vector<int>& mask,
                  const std::vector<POD>& expected_global, gidx_t expected_dof) {
    EXPECT_EQ(gather_scatter.glb_dof(), expected_dof);

    for (int root = 0; root < 3; ++root) {
        std::vector<POD> local(global_index.size(), -1.);
        for (idx_t index = 0; index < local.size(); ++index) {
            if (!mask[index]) {
                local[index] = static_cast<POD>(global_index[index] * 10);
            }
        }

        std::vector<POD> global(fixture.rank == root ? expected_global.size() : 0, -1.);
        idx_t strides[] = {1};
        idx_t extents[] = {1};
        gather_scatter.gather(local.data(), strides, extents, 1, global.data(), strides, extents, 1, root);
        if (fixture.rank == root) {
            EXPECT(global == expected_global);
        }

        local.assign(local.size(), -1.);
        gather_scatter.scatter(global.data(), strides, extents, 1, local.data(), strides, extents, 1, root);
        for (idx_t index = 0; index < local.size(); ++index) {
            if (!mask[index]) {
                EXPECT_EQ(local[index], static_cast<POD>(global_index[index] * 10));
            }
            else {
                EXPECT_EQ(local[index], -1.);
            }
        }
    }
}

//-----------------------------------------------------------------------------

CASE("test_gather") {
    Fixture f;

    SECTION("test_gather_rank0") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            std::vector<POD> loc(f.Nl);
            std::vector<POD> glb(f.Ng());

            for (int j = 0; j < f.Nl; ++j) {
                loc[j] = (idx_t(f.part[j]) != f.rank ? 0 : f.gidx[j] * 10);
            }

            idx_t strides[] = {1};
            idx_t extents[] = {1};
            f.gather_scatter.gather(loc.data(), strides, extents, 1, glb.data(), strides, extents, 1, f.root);

            if (f.rank == f.root) {
                POD glb_c[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
                EXPECT(glb == eckit::testing::make_view(glb_c, glb_c + f.Ng()));
            }
        }
    }

#if 1
    SECTION("test_gather_rank1_deprecated") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl, 2);
            array::ArrayT<POD> glb(f.Ng(), 2);
            array::ArrayT<POD> glb1(f.Ng(), 1);
            array::ArrayT<POD> glb2(f.Ng(), 1);
            array::ArrayView<POD, 2> locv = array::make_view<POD, 2>(loc);
            for (int j = 0; j < f.Nl; ++j) {
                locv(j, 0) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
                locv(j, 1) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
            }

            // Gather complete field
            {
                idx_t loc_strides[] = {loc.stride(0), loc.stride(1)};
                idx_t loc_extents[] = {1, loc.shape(1)};
                idx_t glb_strides[] = {glb.stride(0), glb.stride(1)};
                idx_t glb_extents[] = {1, glb.shape(1)};

                f.gather_scatter.gather(loc.data<POD>(), loc_strides, loc_extents, 2, glb.data<POD>(), glb_strides,
                                        glb_extents, 2, f.root);
            }
            if (f.rank == f.root) {
                auto glbv   = array::make_view<POD, 2>(glb);
                POD glb_c[] = {10, 100, 20, 200, 30, 300, 40, 400, 50, 500, 60, 600, 70, 700, 80, 800, 90, 900};
                idx_t c(0);
                for (idx_t i = 0; i < glb.shape(0); ++i) {
                    for (idx_t j = 0; j < glb.shape(1); ++j) {
                        EXPECT(glbv(i, j) == glb_c[c++]);
                    }
                }
            }

            // Gather only first component
            {
                idx_t loc_strides[] = {loc.stride(0), 2};
                idx_t loc_extents[] = {1, 1};
                idx_t glb_strides[] = {glb1.stride(0), glb1.stride(1)};
                idx_t glb_extents[] = {1, glb1.shape(1)};

                f.gather_scatter.gather(loc.data<POD>(), loc_strides, loc_extents, 2, glb1.data<POD>(), glb_strides,
                                        glb_extents, 2, f.root);
            }
            if (f.rank == f.root) {
                auto glbv    = array::make_view<POD, 2>(glb1);
                POD glb1_c[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
                idx_t c(0);
                for (idx_t i = 0; i < glb1.shape(0); ++i) {
                    for (idx_t j = 0; j < glb1.shape(1); ++j) {
                        EXPECT(glbv(i, j) == glb1_c[c++]);
                    }
                }
            }

            // Gather only second component
            {
                idx_t loc_strides[] = {loc.stride(0), 2};
                idx_t loc_extents[] = {1, 1};
                idx_t glb_strides[] = {glb2.stride(0), glb2.stride(1)};
                idx_t glb_extents[] = {1, glb2.shape(1)};
                f.gather_scatter.gather(loc.data<POD>() + 1, loc_strides, loc_extents, 1, glb2.data<POD>(), glb_strides,
                                        glb_extents, 1, f.root);
            }
            if (f.rank == f.root) {
                auto glbv    = array::make_view<POD, 2>(glb2);
                POD glb2_c[] = {100, 200, 300, 400, 500, 600, 700, 800, 900};
                idx_t c(0);
                for (idx_t i = 0; i < glb2.shape(0); ++i) {
                    for (idx_t j = 0; j < glb2.shape(1); ++j) {
                        EXPECT(glbv(i, j) == glb2_c[c++]);
                    }
                }
            }
        }
    }
#endif

    SECTION("test_gather_rank1") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl, 2);
            array::ArrayT<POD> glb(f.Ng(), 2);
            array::ArrayT<POD> glb1(f.Ng(), 1);
            array::ArrayT<POD> glb2(f.Ng(), 1);
            array::ArrayView<POD, 2> locv = array::make_view<POD, 2>(loc);
            for (int j = 0; j < f.Nl; ++j) {
                locv(j, 0) = (idx_t(f.part[j]) != f.rank ? 0 : f.gidx[j] * 10);
                locv(j, 1) = (idx_t(f.part[j]) != f.rank ? 0 : f.gidx[j] * 100);
            }

// Gather complete field
#if 0
        {
        idx_t loc_strides[] = {2,1};
        idx_t loc_extents[] = {idx_t(f.Nl), 2};
        idx_t loc_rank = 2;
        idx_t loc_mpl_idxpos[] = {0};
        idx_t loc_mpl_rank = 1;
        idx_t glb_strides[] = {2,1};
        idx_t glb_extents[] = {idx_t(f.Ng()), 2};
        idx_t glb_rank = 2;
        idx_t glb_mpl_idxpos[] = {0};
        idx_t glb_mpl_rank = 1;
        parallel::detail::MPL_ArrayView<POD> lview(loc.data<POD>(),loc_strides,loc_extents,loc_rank,loc_mpl_idxpos,loc_mpl_rank);
        parallel::detail::MPL_ArrayView<POD> gview(glb.data(),glb_strides,glb_extents,glb_rank,glb_mpl_idxpos,glb_mpl_rank);

        EXPECT(lview.var_rank() == 1);
        EXPECT(lview.var_stride(0) == 1);
        EXPECT(lview.var_shape(0) == 2);
        EXPECT(gview.var_rank() == 1);
        EXPECT(gview.var_stride(0) == 1);
        EXPECT(gview.var_shape(0) == 2);

        EXPECT(lview.mpl_rank() == 1);
        EXPECT(lview.mpl_stride(0) == 2);
        EXPECT(lview.mpl_shape(0) == f.Nl);
        EXPECT(gview.mpl_rank() == 1);
        EXPECT(gview.mpl_stride(0) == 2);
        EXPECT(gview.mpl_shape(0) == f.Ng());

        f.gather_scatter.gather( loc.data<POD>(), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                              glb.data<POD>(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                              f.root );

        if( mpi::comm().rank() == f.root )
          {
            POD glb_c[] = { 10,100, 20,200, 30,300, 40,400, 50,500, 60,600, 70,700, 80,800, 90,900 };
            EXPECT(make_view(glb.data<POD>(),glb.data<POD>()+2*f.Ng()) == make_view(glb_c,glb_c+2*f.Ng()));
          }
        }
#endif

// Gather only first component
#if 0
        {
          idx_t loc_strides[] = {2,2};
          idx_t loc_extents[] = {idx_t(f.Nl), 1};
          idx_t loc_rank = 2;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {1};
          idx_t glb_extents[] = {idx_t(f.Ng())};
          idx_t glb_rank = 1;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;
          parallel::detail::MPL_ArrayView<POD> lview(loc.data<POD>(),loc_strides,loc_extents,loc_rank,loc_mpl_idxpos,loc_mpl_rank);
          EXPECT(lview.var_rank() == 1);
          EXPECT(lview.var_stride(0) == 2);
          EXPECT(lview.var_shape(0) == 1);
          EXPECT(lview.mpl_rank() == 1);
          EXPECT(lview.mpl_stride(0) == 2);
          EXPECT(lview.mpl_shape(0) == f.Nl);

          parallel::detail::MPL_ArrayView<POD> gview(glb1.data<POD>(),glb_strides,glb_extents,glb_rank,glb_mpl_idxpos,glb_mpl_rank);
          EXPECT(gview.var_rank() == 1);
          EXPECT(gview.var_stride(0) == 1);
          EXPECT(gview.var_shape(0) == 1);
          EXPECT(gview.mpl_rank() == 1);
          EXPECT(gview.mpl_stride(0) == 1);
          EXPECT(gview.mpl_shape(0) == f.Ng());

          f.gather_scatter.gather( loc.data<POD>(),  loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glb1.data<POD>(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        if( mpi::comm().rank() == f.root )
          {
            POD glb1_c[] = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };
            EXPECT(make_view(glb1.data<POD>(),glb1.data<POD>()+f.Ng()) == make_view( glb1_c,glb1_c+f.Ng()));
          }
        }
#endif

// Gather only second component
#if 0
        {
          idx_t loc_strides[] = {2,2};
          idx_t loc_extents[] = {idxt(f.Nl), 1};
          idx_t loc_rank = 2;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {1};
          idx_t glb_extents[] = {idx_t(f.Ng())};
          idx_t glb_rank = 1;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;
          f.gather_scatter.gather( loc.data<POD>()+1,  loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glb2.data<POD>(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb2_c[] = { 100, 200, 300, 400, 500, 600, 700, 800, 900 };
          EXPECT(make_view(glb2.data<POD>(),glb2.data<POD>()+f.Ng()) == make_view(glb2_c,glb2_c+f.Ng()));
        }
#endif
        }
    }

    SECTION("test_gather_rank2") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl, 3, 2);
            array::ArrayT<POD> glb(f.Ng(), 3, 2);
            array::ArrayT<POD> glbx1(f.Ng(), 3);
            array::ArrayT<POD> glbx2(f.Ng(), 3);
            array::ArrayT<POD> glb1x(f.Ng(), 2);
            array::ArrayT<POD> glb2x(f.Ng(), 2);
            array::ArrayT<POD> glb32(f.Ng());

            array::ArrayView<POD, 3> locv = array::make_view<POD, 3>(loc);
            for (int p = 0; p < f.Nl; ++p) {
                for (int i = 0; i < 3; ++i) {
                    locv(p, i, 0) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
                    locv(p, i, 1) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
                }
            }

// Gather complete field
#if 0
        {
          idx_t loc_strides[] = {6,2,1};
          idx_t loc_extents[] = {idx_t(f.Nl), 3, 2};
          idx_t loc_rank = 3;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {6,2,1};
          idx_t glb_extents[] = {idx_t(f.Ng()), 3, 2};
          idx_t glb_rank = 3;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;
          f.gather_scatter.gather( loc.data<POD>(), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glb.data<POD>(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb_c[] = { -1,1, -10,10, -100,100,
                          -2,2, -20,20, -200,200,
                          -3,3, -30,30, -300,300,
                          -4,4, -40,40, -400,400,
                          -5,5, -50,50, -500,500,
                          -6,6, -60,60, -600,600,
                          -7,7, -70,70, -700,700,
                          -8,8, -80,80, -800,800,
                          -9,9, -90,90, -900,900 };
          EXPECT(make_view(glb.data(),glb.data()+6*f.Ng()) == make_view(glb_c,glb_c+6*f.Ng()));
        }
#endif

// Gather var 1
#if 0
        {
          idx_t loc_strides[] = {6,2,2};
          idx_t loc_extents[] = {idx_t(f.Nl), 3, 1};
          idx_t loc_rank = 3;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {6,1};
          idx_t glb_extents[] = {idx_t(f.Ng()), 3};
          idx_t glb_rank = 2;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;
          f.gather_scatter.gather( &locv(0,0,0), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glbx1.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb_c[] = { -1, -10, -100,
                          -2, -20, -200,
                          -3, -30, -300,
                          -4, -40, -400,
                          -5, -50, -500,
                          -6, -60, -600,
                          -7, -70, -700,
                          -8, -80, -800,
                          -9, -90, -900 };
          EXPECT(make_view(glbx1.data(),glbx1.data()+3*f.Ng()) == make_view(glb_c,glb_c+3*f.Ng()));
        }
#endif

// Gather var 2
#if 0
        {
          idx_t loc_strides[] = {6,2,2};
          idx_t loc_extents[] = {idx_t(f.Nl), 3, 1};
          idx_t loc_rank = 3;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {6,1};
          idx_t glb_extents[] = {idx_t(f.Ng()), 3};
          idx_t glb_rank = 2;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;
          f.gather_scatter.gather( &locv(0,0,1), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glbx2.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb_c[] = { 1, 10, 100,
                          2, 20, 200,
                          3, 30, 300,
                          4, 40, 400,
                          5, 50, 500,
                          6, 60, 600,
                          7, 70, 700,
                          8, 80, 800,
                          9, 90, 900 };
          BOOST_CHECK_EQUAL_COLLECTIONS(glbx2.data(),glbx2.data()+3*f.Ng(), glb_c,glb_c+3*f.Ng());
        }
#endif

// Gather lev 1
#if 0
        {
          idx_t loc_strides[] = {6,6,1};
          idx_t loc_extents[] = {idx_t(f.Nl), 1, 2};
          idx_t loc_rank = 3;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {2,1};
          idx_t glb_extents[] = {idx_t(f.Ng()), 2};
          idx_t glb_rank = 2;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;

          f.gather_scatter.gather( &locv(0,0,0), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glb1x.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb_c[] = { -1,1,
                          -2,2,
                          -3,3,
                          -4,4,
                          -5,5,
                          -6,6,
                          -7,7,
                          -8,8,
                          -9,9, };
          BOOST_CHECK_EQUAL_COLLECTIONS(glb1x.data(),glb1x.data()+2*f.Ng(), glb_c,glb_c+2*f.Ng());
        }
#endif

// Gather lev 2
#if 0
        {
          idx_t loc_strides[] = {6,6,1};
          idx_t loc_extents[] = {idx_t(f.Nl), 1, 2};
          idx_t loc_rank = 3;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {2,1};
          idx_t glb_extents[] = {idx_t(f.Ng()), 2};
          idx_t glb_rank = 2;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;
          f.gather_scatter.gather( &locv(0,1,0), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glb2x.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb_c[] = { -10,10,
                          -20,20,
                          -30,30,
                          -40,40,
                          -50,50,
                          -60,60,
                          -70,70,
                          -80,80,
                          -90,90, };
          BOOST_CHECK_EQUAL_COLLECTIONS(glb2x.data(),glb2x.data()+2*f.Ng(), glb_c,glb_c+2*f.Ng());
        }
#endif

// Gather lev 3 var 2
#if 0
        {
          idx_t loc_strides[] = {6,6,2};
          idx_t loc_extents[] = {idx_t(f.Nl), 1, 1};
          idx_t loc_rank = 3;
          idx_t loc_mpl_idxpos[] = {0};
          idx_t loc_mpl_rank = 1;
          idx_t glb_strides[] = {1};
          idx_t glb_extents[] = {idx_t(f.Ng())};
          idx_t glb_rank = 1;
          idx_t glb_mpl_idxpos[] = {0};
          idx_t glb_mpl_rank = 1;

          f.gather_scatter.gather( &locv(0,2,1), loc_strides, loc_extents, loc_rank, loc_mpl_idxpos, loc_mpl_rank,
                                glb32.data(), glb_strides, glb_extents, glb_rank, glb_mpl_idxpos, glb_mpl_rank,
                                f.root );
        }
        if( mpi::comm().rank() == f.root )
        {
          POD glb_c[] = { 100,
                          200,
                          300,
                          400,
                          500,
                          600,
                          700,
                          800,
                          900 };
          BOOST_CHECK_EQUAL_COLLECTIONS(glb32.data(),glb32.data()+f.Ng(), glb_c,glb_c+f.Ng());
        }
#endif
        }
    }

    SECTION("test_gather_rank0_ArrayView") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl);
            array::ArrayT<POD> glb(f.Ng());

            array::ArrayView<POD, 1> locv = array::make_view<POD, 1>(loc);
            array::ArrayView<POD, 1> glbv = array::make_view<POD, 1>(glb);
            for (int p = 0; p < f.Nl; ++p) {
                locv(p) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * 10);
            }

            // Gather complete field
            { f.gather_scatter.gather(locv, glbv, f.root); }
            if (f.rank == f.root) {
                POD glb_c[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
                for (idx_t n = 0; n < glb.shape(0); ++n) {
                    EXPECT(glbv(n) == glb_c[n]);
                }
            }
        }
    }

    SECTION("test_gather_rank1_ArrayView") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl, 2);
            array::ArrayT<POD> glb(f.Ng(), 2);

            array::ArrayView<POD, 2> locv = array::make_view<POD, 2>(loc);
            array::ArrayView<POD, 2> glbv = array::make_view<POD, 2>(glb);
            for (int p = 0; p < f.Nl; ++p) {
                locv(p, 0) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * 10);
                locv(p, 1) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * 10);
            }

            // Gather complete field
            { f.gather_scatter.gather(locv, glbv, f.root); }
            if (f.rank == f.root) {
                POD glb_c[] = {-10, 10, -20, 20, -30, 30, -40, 40, -50, 50, -60, 60, -70, 70, -80, 80, -90, 90};

                auto glbv = array::make_view<POD, 2>(glb);
                idx_t c(0);
                for (idx_t i = 0; i < glb.shape(0); ++i) {
                    for (idx_t j = 0; j < glb.shape(1); ++j) {
                        EXPECT(glbv(i, j) == glb_c[c++]);
                    }
                }
            }
        }
    }

    SECTION("test_gather_rank2_ArrayView") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl, 3, 2);
            array::ArrayT<POD> glb(f.Ng(), 3, 2);

            array::ArrayView<POD, 3> locv = array::make_view<POD, 3>(loc);
            array::ArrayView<POD, 3> glbv = array::make_view<POD, 3>(glb);
            for (int p = 0; p < f.Nl; ++p) {
                for (int i = 0; i < 3; ++i) {
                    locv(p, i, 0) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
                    locv(p, i, 1) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
                }
            }

            // Gather complete field
            { f.gather_scatter.gather(locv, glbv, f.root); }
            if (f.rank == f.root) {
                POD glb_c[] = {-1, 1, -10, 10, -100, 100, -2, 2, -20, 20, -200, 200, -3, 3, -30, 30, -300, 300,
                               -4, 4, -40, 40, -400, 400, -5, 5, -50, 50, -500, 500, -6, 6, -60, 60, -600, 600,
                               -7, 7, -70, 70, -700, 700, -8, 8, -80, 80, -800, 800, -9, 9, -90, 90, -900, 900};
                idx_t c(0);
                for (idx_t i = 0; i < glb.shape(0); ++i) {
                    for (idx_t j = 0; j < glb.shape(1); ++j) {
                        for (idx_t k = 0; k < glb.shape(2); ++k) {
                            EXPECT(glbv(i, j, k) == glb_c[c++]);
                        }
                    }
                }
            }
        }
        f.root = 0;
    }

    SECTION("test_scatter_rank2_ArrayView") {
        for (f.root = 0; f.root < f.comm_size; ++f.root) {
            array::ArrayT<POD> loc(f.Nl, 3, 2);
            array::ArrayT<POD> glb(f.Ng(), 3, 2);

            array::ArrayView<POD, 3> locv = array::make_view<POD, 3>(loc);
            array::ArrayView<POD, 3> glbv = array::make_view<POD, 3>(glb);
            if (f.rank == f.root) {
                POD glb_c[] = {-1, 1, -10, 10, -100, 100, -2, 2, -20, 20, -200, 200, -3, 3, -30, 30, -300, 300,
                               -4, 4, -40, 40, -400, 400, -5, 5, -50, 50, -500, 500, -6, 6, -60, 60, -600, 600,
                               -7, 7, -70, 70, -700, 700, -8, 8, -80, 80, -800, 800, -9, 9, -90, 90, -900, 900};
                idx_t c(0);
                for (int i = 0; i < glb.shape(0); ++i) {
                    for (int j = 0; j < glb.shape(1); ++j) {
                        for (int k = 0; k < glb.shape(2); ++k) {
                            glbv(i, j, k) = glb_c[c++];
                        }
                    }
                }
            }

            POD nan = -1000.;
            locv.assign(nan);

            f.gather_scatter.scatter(glbv, locv, f.root);

            switch (mpi::comm().rank()) {
                case 0: {
                    POD loc_c[] = {nan, nan, nan, nan, nan,  nan, -1,  1,   -10, 10,  -100, 100,
                                   -2,  2,   -20, 20,  -200, 200, -3,  3,   -30, 30,  -300, 300,
                                   nan, nan, nan, nan, nan,  nan, nan, nan, nan, nan, nan,  nan};

                    idx_t c(0);
                    for (idx_t i = 0; i < loc.shape(0); ++i) {
                        for (idx_t j = 0; j < loc.shape(1); ++j) {
                            for (idx_t k = 0; k < loc.shape(2); ++k) {
                                EXPECT(is_approximately_equal(locv(i, j, k), loc_c[c++]));
                            }
                        }
                    }
                    break;
                }
                case 1: {
                    POD loc_c[] = {nan, nan, nan, nan, nan,  nan, -4,  4,   -40, 40,  -400, 400,
                                   -5,  5,   -50, 50,  -500, 500, -6,  6,   -60, 60,  -600, 600,
                                   nan, nan, nan, nan, nan,  nan, nan, nan, nan, nan, nan,  nan};
                    idx_t c(0);
                    for (idx_t i = 0; i < loc.shape(0); ++i) {
                        for (idx_t j = 0; j < loc.shape(1); ++j) {
                            for (idx_t k = 0; k < loc.shape(2); ++k) {
                                EXPECT(is_approximately_equal(locv(i, j, k), loc_c[c++]));
                            }
                        }
                    }
                    break;
                }
                case 2: {
                    POD loc_c[] = {nan,  nan, nan,  nan, nan, nan, nan, nan, nan,  nan, nan, nan, -7,  7,
                                   -70,  70,  -700, 700, -8,  8,   -80, 80,  -800, 800, -9,  9,   -90, 90,
                                   -900, 900, nan,  nan, nan, nan, nan, nan, nan,  nan, nan, nan, nan, nan};
                    idx_t c(0);
                    for (idx_t i = 0; i < loc.shape(0); ++i) {
                        for (idx_t j = 0; j < loc.shape(1); ++j) {
                            for (idx_t k = 0; k < loc.shape(2); ++k) {
                                EXPECT(is_approximately_equal(locv(i, j, k), loc_c[c++]));
                            }
                        }
                    }
                    break;
                }
            }
        }
    }
}

CASE("test_gather_sparse_global_indices") {
    SparseGlobalIndexFixture f;
    // The fixture includes global-index 14399999990.
    // An 64-bit-value array of this size (a global gathered field) would
    // require more than 100GB.
    constexpr size_t sparse_global_index_memory_limit = 10 * GB; // less than 100 GB

    for (f.root = 0; f.root < f.comm_size; ++f.root) {
        std::vector<POD> loc(f.Nl);
        std::vector<POD> glb(f.Ng());

        for (int j = 0; j < f.Nl; ++j) {
            loc[j] = (idx_t(f.part[j]) != f.rank ? -1. : f.gidx[j] * 10.);
        }

        idx_t strides[] = {1};
        idx_t extents[] = {1};
        f.gather_scatter.gather(loc.data(), strides, extents, 1, glb.data(), strides, extents, 1, f.root);

        EXPECT(f.gather_scatter.glb_dof() == 6);

        if (f.rank == f.root) {
            POD glb_c[] = {100, 200, 300, 400, 500, 14399999990.};
            EXPECT(glb == eckit::testing::make_view(glb_c, glb_c + f.Ng()));
        }
    }

    Log::warning() << "test_gather_sparse_global_indices: peak memory usage = " << peakMemory() / MB << " MB" << std::endl;
    EXPECT(peakMemory() < sparse_global_index_memory_limit);
}

CASE("test_gather_setup_options") {
    SetupOptionsFixture fixture;
    const idx_t size = static_cast<idx_t>(fixture.part.size());
    const std::vector<POD> dense_expected{10., 20., 30., 40., 50., 60.};
    const std::vector<POD> sparse_compacted_expected{20., 40., 70., 90., 110., 130.};
    const std::vector<POD> sparse_preserved_expected{-1., 20., -1., 40., -1., -1., 70.,
                                                      -1., 90., -1., 110., -1., 130.};

    SECTION("automatic mask with default communicator") {
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(fixture.part.data(), fixture.remote_idx.data(), 0, fixture.dense_gidx.data(), size);
        test_mapping(gather_scatter, fixture, fixture.dense_gidx, fixture.mask, dense_expected, 6);
    }

    SECTION("automatic mask with named communicator") {
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), fixture.remote_idx.data(), 0,
                             fixture.dense_gidx.data(), size);
        test_mapping(gather_scatter, fixture, fixture.dense_gidx, fixture.mask, dense_expected, 6);
    }

    SECTION("explicit mask with default communicator") {
        const std::vector<int> mask{0, 1, 1, 1};
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(fixture.part.data(), fixture.remote_idx.data(), 0, fixture.dense_gidx.data(),
                             mask.data(), size);
        test_mapping(gather_scatter, fixture, fixture.dense_gidx, mask, {10., 20., 30.}, 3);
    }

    SECTION("explicit mask with named communicator") {
        const std::vector<int> mask{0, 1, 1, 1};
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), fixture.remote_idx.data(), 0,
                             fixture.dense_gidx.data(), mask.data(), size);
        test_mapping(gather_scatter, fixture, fixture.dense_gidx, mask, {10., 20., 30.}, 3);
    }

    SECTION("remote index base one") {
        auto remote_idx = fixture.remote_idx;
        for (auto& index : remote_idx) {
            ++index;
        }
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), remote_idx.data(), 1,
                             fixture.dense_gidx.data(), fixture.mask.data(), size);
        test_mapping(gather_scatter, fixture, fixture.dense_gidx, fixture.mask, dense_expected, 6);
    }

    SECTION("sparse global indices are compacted by default") {
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), fixture.remote_idx.data(), 0,
                             fixture.sparse_gidx.data(), fixture.mask.data(), size, false);
        test_mapping(gather_scatter, fixture, fixture.sparse_gidx, fixture.mask, sparse_compacted_expected, 6);
    }

    SECTION("sparse global indices can be preserved") {
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), fixture.remote_idx.data(), 0,
                             fixture.sparse_gidx.data(), fixture.mask.data(), size, true);
        test_mapping(gather_scatter, fixture, fixture.sparse_gidx, fixture.mask, sparse_preserved_expected, 6);
    }

    SECTION("duplicate global indices are deduplicated when compacting") {
        const std::vector<int> mask(size, 0);
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), fixture.remote_idx.data(), 0,
                             fixture.dense_gidx.data(), mask.data(), size, false);
        EXPECT_EQ(gather_scatter.glb_dof(), 6);

        for (int root = 0; root < 3; ++root) {
            std::vector<POD> local(size);
            for (idx_t index = 0; index < size; ++index) {
                local[index] = static_cast<POD>(fixture.dense_gidx[index] * 10);
            }
            std::vector<POD> global(fixture.rank == root ? dense_expected.size() : 0, -1.);
            idx_t strides[] = {1};
            idx_t extents[] = {1};
            gather_scatter.gather(local.data(), strides, extents, 1, global.data(), strides, extents, 1, root);
            if (fixture.rank == root) {
                EXPECT(global == dense_expected);
            }
        }
    }

    SECTION("all points masked") {
        std::vector<int> mask(size, 1);
        parallel::GatherScatter gather_scatter;
        gather_scatter.setup(mpi::comm().name(), fixture.part.data(), fixture.remote_idx.data(), 0,
                             fixture.dense_gidx.data(), mask.data(), size);
        EXPECT_EQ(gather_scatter.glb_dof(), 0);
        EXPECT_EQ(gather_scatter.loc_dof(), 0);
    }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
