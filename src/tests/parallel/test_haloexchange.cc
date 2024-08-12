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
#include <memory>
#include <sstream>

#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/MakeView.h"
#include "atlas/library/config.h"
#include "atlas/parallel/HaloExchange.h"
#include "atlas/parallel/mpi/mpi.h"

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

template <int Rank, typename FirstDim>
size_t eval_idx(size_t pos, std::array<size_t, Rank>& strides, FirstDim first) {
    return first * strides[pos];
}

template <int Rank, typename FirstDim, typename SecondDim>
size_t eval_idx(size_t pos, std::array<size_t, Rank>& strides, FirstDim first, SecondDim second) {
    return first * strides[pos] + eval_idx<Rank>(pos + 1, strides, second);
}

template <int Rank, typename FirstDim, typename SecondDim, typename ThirdDim>
size_t eval_idx(size_t pos, std::array<size_t, Rank>& strides, FirstDim first, SecondDim second, ThirdDim third) {
    return first * strides[pos] + eval_idx<Rank>(pos + 1, strides, second, third);
}

template <typename DATA_TYPE, int Rank, int Dim>
struct validate_impl;

template <typename DATA_TYPE, int Rank>
struct validate_impl<DATA_TYPE, Rank, 0> {
    template <typename... Int>
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv, DATA_TYPE arr_c[], std::array<size_t, Rank>& strides,
                      Int... dims) {
        EXPECT(arrv(dims...) == arr_c[eval_idx<Rank>((size_t)0, strides, dims...)]);
    }
};

template <typename DATA_TYPE, int Rank, int Dim>
struct validate_impl {
    template <typename... Int>
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv, DATA_TYPE arr_c[], std::array<size_t, Rank>& strides,
                      Int... dims) {
        for (idx_t cnt = 0; cnt < arrv.template shape<Rank - Dim>(); ++cnt) {
            validate_impl<DATA_TYPE, Rank, Dim - 1>::apply(arrv, arr_c, strides, dims..., cnt);
        }
    }
};

template <int Dim>
struct compute_strides;

template <>
struct compute_strides<0> {
    template <typename DATA_TYPE, int Rank>
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv, std::array<size_t, (unsigned int)Rank>& strides) {}
};

template <int Dim>
struct compute_strides {
    template <typename DATA_TYPE, int Rank>
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv, std::array<size_t, (unsigned int)Rank>& strides) {
        strides[Dim - 1] = strides[Dim] * arrv.template shape<(unsigned int)Dim>();
        compute_strides<Dim - 1>::apply(arrv, strides);
    }
};

template <typename DATA_TYPE, int Rank>
struct validate {
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv, DATA_TYPE arr_c[]) {
        std::array<size_t, Rank> strides;
        strides[Rank - 1] = 1;
        compute_strides<Rank - 1>::apply(arrv, strides);

        for (idx_t i = 0; i < arrv.template shape<0>(); ++i) {
            validate_impl<DATA_TYPE, Rank, Rank - 1>::apply(arrv, arr_c, strides, i);
        }
    }
};

struct Fixture {
    Fixture(bool on_device = false): on_device_(on_device) {
        int nnodes_c[] = {5, 6, 7};
        nb_nodes       = vec(nnodes_c);
        N              = nb_nodes[mpi::comm().rank()];
        switch (mpi::comm().rank()) {
            case 0: {
                int part_c[]   = {2, 0, 0, 0, 1};
                part           = vec(part_c);
                idx_t ridx_c[] = {4, 1, 2, 3, 1};
                ridx           = vec(ridx_c);
                POD gidx_c[]   = {0, 1, 2, 3, 0};
                gidx           = vec(gidx_c);
                break;
            }
            case 1: {
                int part_c[]   = {0, 1, 1, 1, 2, 2};
                part           = vec(part_c);
                idx_t ridx_c[] = {3, 1, 2, 3, 2, 3};
                ridx           = vec(ridx_c);
                POD gidx_c[]   = {0, 4, 5, 6, 0, 0};
                gidx           = vec(gidx_c);
                break;
            }
            case 2: {
                int part_c[]   = {1, 1, 2, 2, 2, 0, 0};
                part           = vec(part_c);
                idx_t ridx_c[] = {2, 3, 2, 3, 4, 1, 2};
                ridx           = vec(ridx_c);
                POD gidx_c[]   = {0, 0, 7, 8, 9, 0, 0};
                gidx           = vec(gidx_c);
                break;
            }
        }
        halo_exchange.setup(part.data(), ridx.data(), 0, N);
    }
    parallel::HaloExchange halo_exchange;
    std::vector<int> nb_nodes;
    std::vector<int> part;
    std::vector<idx_t> ridx;
    std::vector<POD> gidx;

    int N;
    bool on_device_;
};

//-----------------------------------------------------------------------------

void test_rank0_arrview(Fixture& f) {
    array::ArrayT<POD> arr(f.N);
    array::ArrayView<POD, 1> arrv = array::make_host_view<POD, 1>(arr);
    for (int j = 0; j < f.N; ++j) {
        arrv(j) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j]);
    }

    arr.syncHostDevice();

    f.halo_exchange.execute<POD, 1>(arr, f.on_device_);

    arr.syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank1(Fixture& f) {
    array::ArrayT<POD> arr(f.N, 2);
    array::ArrayView<POD, 2> arrv = array::make_host_view<POD, 2>(arr);
    for (int j = 0; j < f.N; ++j) {
        arrv(j, 0) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv(j, 1) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
    }

    arr.syncHostDevice();

    f.halo_exchange.execute<POD, 2>(arr, f.on_device_);

    arr.syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {90, 900, 10, 100, 20, 200, 30, 300, 40, 400};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {30, 300, 40, 400, 50, 500, 60, 600, 70, 700, 80, 800};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {50, 500, 60, 600, 70, 700, 80, 800, 90, 900, 10, 100, 20, 200};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank1_strided_v1(Fixture& f) {
    // create a 2d field from the gidx data, with two components per grid point
    array::ArrayT<POD> arr_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_t = array::make_host_view<POD, 2>(arr_t);
    for (int j = 0; j < f.N; ++j) {
        arrv_t(j, 0) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv_t(j, 1) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
    }

    // create a wrap array where we fake the strides in a way that the second
    // dimension
    // (number of components) contains only one component but the associated
    // stride is 2
    // (i.e. we are only selecting and exchanging the first component of the
    // field)

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        arrv_t.data(),
        array::ArraySpec {
            array::make_shape(f.N, 1),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
            array::make_strides(32, 1)
#else
            array::make_strides(2, 1)
#endif
        }));

    arr->syncHostDevice();

    f.halo_exchange.execute<POD, 2>(*arr, f.on_device_);

    arr->syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {90, 0, 10, 100, 20, 200, 30, 300, 40, 0};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {30, 0, 40, 400, 50, 500, 60, 600, 70, 0, 80, 0};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {50, 0, 60, 0, 70, 700, 80, 800, 90, 900, 10, 0, 20, 0};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
    }
}

void test_rank1_strided_v2(Fixture& f) {
    // create a 2d field from the gidx data, with two components per grid point
    array::ArrayT<POD> arr_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_t = array::make_host_view<POD, 2>(arr_t);
    for (int j = 0; j < f.N; ++j) {
        arrv_t(j, 0) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv_t(j, 1) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
    }

    // create a wrap array where we fake the strides in a way that the second
    // dimension
    // (number of components) contains only one component but the associated
    // stride is 2
    // (i.e. we are only selecting and exchanging the first component of the
    // field)

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        &(arrv_t(0, 1)), array::ArraySpec {
            array::make_shape(f.N, 1),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
            array::make_strides(32, 1)
#else
            array::make_strides(2, 1)
#endif
        }));

    arr->syncHostDevice();

    f.halo_exchange.execute<POD, 2>(*arr, false);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 900, 10, 100, 20, 200, 30, 300, 0, 400};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 300, 40, 400, 50, 500, 60, 600, 0, 700, 0, 800};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 500, 0, 600, 70, 700, 80, 800, 90, 900, 0, 100, 0, 200};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
    }
}

void test_rank2(Fixture& f) {
    array::ArrayT<POD> arr(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv = array::make_host_view<POD, 3>(arr);
    for (int p = 0; p < f.N; ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv((size_t)p, i, (size_t)0) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv((size_t)p, i, (size_t)1) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
        }
    }

    arr.syncHostDevice();

    f.halo_exchange.execute<POD, 3>(arr, f.on_device_);

    arr.syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {-9, 9,    -90, 90, -900, 900, -1, 1,    -10, 10, -100, 100, -2, 2,    -20,
                           20, -200, 200, -3, 3,    -30, 30, -300, 300, -4, 4,    -40, 40, -400, 400};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {-3, 3, -30, 30, -300, 300, -4, 4, -40, 40, -400, 400, -5, 5, -50, 50, -500, 500,
                           -6, 6, -60, 60, -600, 600, -7, 7, -70, 70, -700, 700, -8, 8, -80, 80, -800, 800};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {-5,   5,   -50,  50,  -500, 500, -6,   6,   -60,  60,  -600, 600, -7,   7,
                           -70,  70,  -700, 700, -8,   8,   -80,  80,  -800, 800, -9,   9,   -90,  90,
                           -900, 900, -1,   1,   -10,  10,  -100, 100, -2,   2,   -20,  20,  -200, 200};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank2_l1(Fixture& f) {
    array::ArrayT<POD> arr_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_t = array::make_host_view<POD, 3>(arr_t);
    for (int p = 0; p < f.N; ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv_t((size_t)p, i, (size_t)0) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv_t((size_t)p, i, (size_t)1) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
        }
    }

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        arrv_t.data(), array::ArraySpec {
            array::make_shape(f.N, 1, 2),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
            array::make_strides(96, 32, 1)
#else
            array::make_strides(6, 2, 1)
#endif
        }));

    arr->syncHostDevice();

    f.halo_exchange.execute<POD, 3>(*arr, false);

    arr->syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {-9, 9, 0,   0,  0,    0,    // halo
                           -1, 1, -10, 10, -100, 100,  // core
                           -2, 2, -20, 20, -200, 200,  // core
                           -3, 3, -30, 30, -300, 300,  // core
                           -4, 4, 0,   0,  0,    0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {-3, 3, 0,   0,  0,    0,    // halo
                           -4, 4, -40, 40, -400, 400,  // core
                           -5, 5, -50, 50, -500, 500,  // core
                           -6, 6, -60, 60, -600, 600,  // core
                           -7, 7, 0,   0,  0,    0,    // halo
                           -8, 8, 0,   0,  0,    0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {-5, 5, 0,   0,  0,    0,    // halo
                           -6, 6, 0,   0,  0,    0,    // halo
                           -7, 7, -70, 70, -700, 700,  // core
                           -8, 8, -80, 80, -800, 800,  // core
                           -9, 9, -90, 90, -900, 900,  // core
                           -1, 1, 0,   0,  0,    0,    // halo
                           -2, 2, 0,   0,  0,    0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
    }
}

void test_rank2_l2_v2(Fixture& f) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    // Test rank 2 halo-exchange
    array::ArrayT<POD> arr_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_t = array::make_host_view<POD, 3>(arr_t);
    for (int p = 0; p < f.N; ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv_t((size_t)p, i, (size_t)0) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv_t((size_t)p, i, (size_t)1) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
        }
    }

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        &arrv_t(0, 1, 1), array::ArraySpec {
            array::make_shape(f.N, 1, 1),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
            array::make_strides(192, 32, 1)
#else
            array::make_strides(6, 2, 1)
#endif
        }));

    f.halo_exchange.execute<POD, 3>(*arr, f.on_device_);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0,  0, 0,   90, 0,    0,    // halo
                           -1, 1, -10, 10, -100, 100,  // core
                           -2, 2, -20, 20, -200, 200,  // core
                           -3, 3, -30, 30, -300, 300,  // core
                           0,  0, 0,   40, 0,    0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0,  0, 0,   30, 0,    0,    // halo
                           -4, 4, -40, 40, -400, 400,  // core
                           -5, 5, -50, 50, -500, 500,  // core
                           -6, 6, -60, 60, -600, 600,  // core
                           0,  0, 0,   70, 0,    0,    // halo
                           0,  0, 0,   80, 0,    0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0,  0, 0,   50, 0,    0,    // halo
                           0,  0, 0,   60, 0,    0,    // halo
                           -7, 7, -70, 70, -700, 700,  // core
                           -8, 8, -80, 80, -800, 800,  // core
                           -9, 9, -90, 90, -900, 900,  // core
                           0,  0, 0,   10, 0,    0,    // halo
                           0,  0, 0,   20, 0,    0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
    }
#endif
}

void test_rank2_v2(Fixture& f) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    array::ArrayT<POD> arr_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_t = array::make_view<POD, 3>(arr_t);
    for (int p = 0; p < f.N; ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv_t((size_t)p, i, (size_t)0) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv_t((size_t)p, i, (size_t)1) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
        }
    }

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        &arrv_t(0, 0, 1), array::ArraySpec {
            array::make_shape(f.N, 3, 1),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
            array::make_strides(192, 32, 2)
#else
            array::make_strides(6, 2, 2)
#endif
        }));

    f.halo_exchange.execute<POD, 3>(*arr, f.on_device_);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0,  9, 0,   90, 0,    900,   // halo
                           -1, 1, -10, 10, -100, 100,   // core
                           -2, 2, -20, 20, -200, 200,   // core
                           -3, 3, -30, 30, -300, 300,   // core
                           0,  4, 0,   40, 0,    400};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0,  3, 0,   30, 0,    300,   // halo
                           -4, 4, -40, 40, -400, 400,   // core
                           -5, 5, -50, 50, -500, 500,   // core
                           -6, 6, -60, 60, -600, 600,   // core
                           0,  7, 0,   70, 0,    700,   // halo
                           0,  8, 0,   80, 0,    800};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0,  5, 0,   50, 0,    500,   // halo
                           0,  6, 0,   60, 0,    600,   // halo
                           -7, 7, -70, 70, -700, 700,   // core
                           -8, 8, -80, 80, -800, 800,   // core
                           -9, 9, -90, 90, -900, 900,   // core
                           0,  1, 0,   10, 0,    100,   // halo
                           0,  2, 0,   20, 0,    200};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
    }
#endif
}

void test_rank0_wrap(Fixture& f) {
    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(f.gidx.data(), array::make_shape(f.N)));
    array::ArrayView<POD, 1> arrv = array::make_view<POD, 1>(*arr);

    arr->syncHostDevice();

    f.halo_exchange.execute<POD, 1>(*arr, f.on_device_);

    arr->syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank1_paralleldim1(Fixture& f) {
    array::ArrayT<POD> arr(2, f.N);
    array::ArrayView<POD, 2> arrv = array::make_view<POD, 2>(arr);
    for (int j = 0; j < f.N; ++j) {
        arrv(0, j) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv(1, j) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
    }

    f.halo_exchange.execute<POD, 2, array::LastDim>(arr, false);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {90, 10, 20, 30, 40, 900, 100, 200, 300, 400};  // 90,900, 10,100, 20,200, 30,300, 40,400 };
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {30, 40, 50, 60, 70, 80, 300, 400, 500, 600, 700, 800};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {50, 60, 70, 80, 90, 10, 20, 500, 600, 700, 800, 900, 100, 200};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank2_paralleldim2(Fixture& f) {
    array::ArrayT<POD> arr(3, f.N, 2);
    array::ArrayView<POD, 3> arrv = array::make_view<POD, 3>(arr);
    for (int p = 0; p < f.N; ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv(i, (size_t)p, (size_t)0) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv(i, (size_t)p, (size_t)1) = (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
        }
    }

    f.halo_exchange.execute<POD, 3, array::Dim<1>>(arr, false);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {-9, 9,   -1, 1,   -2, 2,    -3,  3,    -4,  4,    -90, 90,   -10, 10,   -20,
                           20, -30, 30, -40, 40, -900, 900, -100, 100, -200, 200, -300, 300, -400, 400};

            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {-3,  3,  -4,  4,  -5,  5,  -6,   6,   -7,   7,   -8,   8,   -30,  30,  -40,  40,  -50,  50,
                           -60, 60, -70, 70, -80, 80, -300, 300, -400, 400, -500, 500, -600, 600, -700, 700, -800, 800};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {-5,   5,   -6,   6,   -7,   7,   -8,   8,   -9,   9,   -1,   1,   -2,   2,
                           -50,  50,  -60,  60,  -70,  70,  -80,  80,  -90,  90,  -10,  10,  -20,  20,
                           -500, 500, -600, 600, -700, 700, -800, 800, -900, 900, -100, 100, -200, 200};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank1_cinterface(Fixture& f) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    array::ArrayT<POD> arr(f.N, 2);
    array::ArrayView<POD, 2> arrv = array::make_host_view<POD, 2>(arr);
    for (int j = 0; j < f.N; ++j) {
        arrv(j, 0) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv(j, 1) = (size_t(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
    }

    arr.syncHostDevice();

    int shapes[2]  = {(int)arrv.shape(0), (int)arrv.shape(1)};
    int strides[2] = {(int)arrv.stride(0), (int)arrv.stride(1)};

    atlas__HaloExchange__execute_strided_double(&(f.halo_exchange), arrv.data(), &(strides[1]), &(shapes[1]), 1);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {90, 900, 10, 100, 20, 200, 30, 300, 40, 400};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {30, 300, 40, 400, 50, 500, 60, 600, 70, 700, 80, 800};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {50, 500, 60, 600, 70, 700, 80, 800, 90, 900, 10, 100, 20, 200};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
    }
#endif
}

CASE("test_haloexchange") {
    Fixture f;

    SECTION("test_rank0_arrview") { test_rank0_arrview(f); }

    SECTION("test_rank1") { test_rank1(f); }

    SECTION("test_rank1_strided_v1") { test_rank1_strided_v1(f); }

    SECTION("test_rank1_strided_v2") { test_rank1_strided_v2(f); }

    SECTION("test_rank2") { test_rank2(f); }

    SECTION("test_rank2_l1") { test_rank2_l1(f); }

    SECTION("test_rank2_l2_v2") { test_rank2_l2_v2(f); }

    SECTION("test_rank2_v2") { test_rank2_v2(f); }

    SECTION("test_rank0_wrap") { test_rank0_wrap(f); }

    SECTION("test_rank1_paralleldim_1") { test_rank1_paralleldim1(f); }

    SECTION("test_rank2_paralleldim_2") { test_rank2_paralleldim2(f); }

    SECTION("test_rank1_cinterface") { test_rank1_cinterface(f); }
}

#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
CASE("test_haloexchange on device") {
    bool on_device = true;
    Fixture f(on_device);

    SECTION("test_rank0_arrview") { test_rank0_arrview(f); }

    SECTION("test_rank1") { test_rank1(f); }

    SECTION("test_rank2") { test_rank2(f); }

    SECTION("test_rank0_wrap") { test_rank0_wrap(f); }
}
#endif


//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
