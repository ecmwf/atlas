/*
 * (C) British Crown Copyright, 2020, Met Office
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

template <typename T, size_t N>
std::vector<T> vec(const T (&list)[N]) {
    return std::vector<T>(list, list + N);
}

template <int Rank, typename FirstDim>
size_t eval_idx(size_t pos, std::array<size_t, Rank>& strides, FirstDim first) {
    return static_cast<size_t>(first) * strides[pos];
}

template <int Rank, typename FirstDim, typename SecondDim>
size_t eval_idx(size_t pos, std::array<size_t, Rank>& strides, FirstDim first, SecondDim second) {
    return static_cast<size_t>(first) * strides[pos] + eval_idx<Rank>(pos + 1, strides, second);
}

template <int Rank, typename FirstDim, typename SecondDim, typename ThirdDim>
size_t eval_idx(size_t pos, std::array<size_t, Rank>& strides, FirstDim first, SecondDim second, ThirdDim third) {
    return static_cast<size_t>(first) * strides[pos] + eval_idx<Rank>(pos + 1, strides, second, third);
}

template <typename DATA_TYPE, int Rank, int Dim>
struct validate_impl;

template <typename DATA_TYPE, int Rank>
struct validate_impl<DATA_TYPE, Rank, 0> {
    template <typename... Int>
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv, DATA_TYPE arr_c[], std::array<size_t, Rank>& strides,
                      Int... dims) {
        EXPECT(arrv(dims...) == arr_c[eval_idx<Rank>(static_cast<size_t>(0), strides, dims...)]);
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
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv,
                      std::array<size_t, static_cast<unsigned int>(Rank)>& strides) {}
};

template <int Dim>
struct compute_strides {
    template <typename DATA_TYPE, int Rank>
    static void apply(array::ArrayView<DATA_TYPE, Rank>& arrv,
                      std::array<size_t, static_cast<unsigned int>(Rank)>& strides) {
        strides[Dim - 1] =
            strides[Dim] * static_cast<unsigned int>(arrv.template shape<static_cast<unsigned int>(Dim)>());
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
    int N;
    std::vector<int> nb_nodes;
    std::vector<int> part;
    std::vector<idx_t> ridx;
    std::vector<POD> gidx;
    bool on_device_;

    std::unique_ptr<parallel::HaloExchange> halo_exchange_std{new parallel::HaloExchange()};

    Fixture(bool on_device): on_device_(on_device) {
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
        halo_exchange_std->setup(part.data(), ridx.data(), 0, N);
    }
};

//-----------------------------------------------------------------------------

void test_rank0_arrview(Fixture& f) {
    array::ArrayT<POD> arr(f.N);
    array::ArrayView<POD, 1> arrv = array::make_host_view<POD, 1>(arr);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {90, 1, 2, 3, 40};
            for (int j = 0; j < f.N; ++j) {
                arrv(j) = arr_c[j];
            }
            break;
        }
        case 1: {
            POD arr_c[] = {30, 4, 5, 6, 70, 80};
            for (int j = 0; j < f.N; ++j) {
                arrv(j) = arr_c[j];
            }
            break;
        }
        case 2: {
            POD arr_c[] = {50, 60, 7, 8, 9, 10, 20};
            for (int j = 0; j < f.N; ++j) {
                arrv(j) = arr_c[j];
            }
            break;
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 1>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 11, 22, 33, 0};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 44, 55, 66, 0, 0};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 77, 88, 99, 0, 0};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
    }
}

//-----------------------------------------------------------------------------

void adjoint_test(POD sum1, POD sum2, const std::string string_test) {
    POD tol(1e-8);

    atlas::mpi::comm().allReduceInPlace(sum1, eckit::mpi::sum());

    atlas::mpi::comm().allReduceInPlace(sum2, eckit::mpi::sum());


    std::cout << "Adjoint test " << string_test << " " << sum1 << " " << sum2 << std::endl;

    EXPECT(std::abs(sum1 - sum2) < tol);
}

//-----------------------------------------------------------------------------
// The *_adj_tests perform the adjoint test of the form
//    < A x , A x> = < A^T A x , x >
//    where A is the halo exchange and A^T is the adjoint of the halo exchange
//    and  < > are the inner products
//-----------------------------------------------------------------------------
void test_rank0_arrview_adj_test(Fixture& f) {
    array::ArrayT<POD> arr_init(f.N);
    array::ArrayT<POD> arr(f.N);
    array::ArrayView<POD, 1> arrv_init = array::make_host_view<POD, 1>(arr_init);
    array::ArrayView<POD, 1> arrv      = array::make_host_view<POD, 1>(arr);

    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        arrv_init(j) = (static_cast<std::size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j]);
        arrv(j)      = arrv_init(j);
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute<POD, 1>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    // sum1
    POD sum1(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        sum1 += arrv(j) * arrv(j);
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 1>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    // sum2
    POD sum2(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        sum2 += arrv_init(j) * arrv(j);
    }

    adjoint_test(sum1, sum2, "test_rank0_arrview_adj_test");
}


void test_rank1(Fixture& f) {
    array::ArrayT<POD> arr(f.N, 2);
    array::ArrayView<POD, 2> arrv = array::make_host_view<POD, 2>(arr);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (int j = 0; j < f.N; ++j) {
                arrv(j, 0) = arr_c[j] * 10;
                arrv(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (int j = 0; j < f.N; ++j) {
                arrv(j, 0) = arr_c[j] * 10;
                arrv(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (int j = 0; j < f.N; ++j) {
                arrv(j, 0) = arr_c[j] * 10;
                arrv(j, 1) = arr_c[j] * 100;
            }
            break;
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 2>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 0, 20, 200, 40, 400, 60, 600, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 0, 80, 800, 100, 1000, 120, 1200, 0, 0, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 0, 0, 140, 1400, 160, 1600, 180, 1800, 0, 0, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank1_adj_test(Fixture& f) {
    array::ArrayT<POD> arr_init(f.N, 2);
    array::ArrayT<POD> arr(f.N, 2);
    array::ArrayView<POD, 2> arrv_init = array::make_host_view<POD, 2>(arr_init);
    array::ArrayView<POD, 2> arrv      = array::make_host_view<POD, 2>(arr);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        arrv_init(j, 0ul) = (static_cast<size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv_init(j, 1ul) = (static_cast<size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
        arrv(j, 0ul)      = arrv_init(j, 0ul);
        arrv(j, 1ul)      = arrv_init(j, 1ul);
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute<POD, 2>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    // sum1
    POD sum1(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum1 += arrv(j, i) * arrv(j, i);
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 2>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    // sum2
    POD sum2(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum2 += arrv_init(j, i) * arrv(j, i);
        }
    }

    adjoint_test(sum1, sum2, "test_rank1_adj_test");
}

void test_rank1_strided_v1(Fixture& f) {
    // create a 2d field from the gidx data, with two components per grid point
    array::ArrayT<POD> arr_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_t = array::make_host_view<POD, 2>(arr_t);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (int j = 0; j < f.N; ++j) {
                arrv_t(j, 0) = arr_c[j] * 10;
                arrv_t(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (int j = 0; j < f.N; ++j) {
                arrv_t(j, 0) = arr_c[j] * 10;
                arrv_t(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (int j = 0; j < f.N; ++j) {
                arrv_t(j, 0) = arr_c[j] * 10;
                arrv_t(j, 1) = arr_c[j] * 100;
            }
            break;
        }
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
        }
#else
                array::make_strides(2, 1)
        }
#endif
        ));

    arr->syncHostDevice();

    f.halo_exchange_std->execute_adjoint<POD, 2>(*arr, f.on_device_);

    arr->syncHostDevice();

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 900, 20, 100, 40, 200, 60, 300, 0, 400};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 300, 80, 400, 100, 500, 120, 600, 0, 700, 0, 800};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 500, 0, 600, 140, 700, 160, 800, 180, 900, 0, 100, 0, 200};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
    }
}

void test_rank1_strided_v1_adj_test(Fixture& f) {
    // create a 2d field from the gidx data, with two components per grid point
    array::ArrayT<POD> arr_init_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_init_t = array::make_host_view<POD, 2>(arr_init_t);
    array::ArrayT<POD> arr_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_t = array::make_host_view<POD, 2>(arr_t);

    for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
        arrv_init_t(j, 0ul) = (static_cast<size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv_init_t(j, 1ul) = (static_cast<size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
        arrv_t(j, 0ul)      = arrv_init_t(j, 0ul);
        arrv_t(j, 1ul)      = arrv_init_t(j, 1ul);
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
        }
#else
                array::make_strides(2, 1)
        }
#endif
        ));

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute<POD, 2>(*arr, f.on_device_);

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    // sum1
    POD sum1(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum1 += arrv_t(j, i) * arrv_t(j, i);
        }
    }

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 2>(*arr, f.on_device_);

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    // sum2
    POD sum2(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum2 += arrv_init_t(j, i) * arrv_t(j, i);
        }
    }

    adjoint_test(sum1, sum2, "test_rank1_strided_v1_adj_test");
}

void test_rank1_strided_v2(Fixture& f) {
    // create a 2d field from the gidx data, with two components per grid point
    array::ArrayT<POD> arr_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_t = array::make_host_view<POD, 2>(arr_t);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (int j = 0; j < f.N; ++j) {
                arrv_t(j, 0) = arr_c[j] * 10;
                arrv_t(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (int j = 0; j < f.N; ++j) {
                arrv_t(j, 0) = arr_c[j] * 10;
                arrv_t(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (int j = 0; j < f.N; ++j) {
                arrv_t(j, 0) = arr_c[j] * 10;
                arrv_t(j, 1) = arr_c[j] * 100;
            }
            break;
        }
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

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 2>(*arr, false);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {90, 0, 10, 200, 20, 400, 30, 600, 40, 0};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {30, 0, 40, 800, 50, 1000, 60, 1200, 70, 0, 80, 0};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {50, 0, 60, 0, 70, 1400, 80, 1600, 90, 1800, 10, 0, 20, 0};
            validate<POD, 2>::apply(arrv_t, arr_c);
            break;
        }
    }
}

void test_rank1_strided_v2_adj_test(Fixture& f) {
    // create a 2d field from the gidx data, with two components per grid point
    array::ArrayT<POD> arr_init_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_init_t = array::make_host_view<POD, 2>(arr_init_t);
    array::ArrayT<POD> arr_t(f.N, 2);
    array::ArrayView<POD, 2> arrv_t = array::make_host_view<POD, 2>(arr_t);
    for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
        arrv_init_t(j, 0ul) = (static_cast<size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv_init_t(j, 1ul) = (static_cast<size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
        arrv_t(j, 0ul)      = arrv_init_t(j, 0ul);
        arrv_t(j, 1ul)      = arrv_init_t(j, 1ul);
    }

    // create a wrap array where we fake the strides in a way that the second
    // dimension
    // (number of components) contains only one component but the associated
    // stride is 2
    // (i.e. we are only selecting and exchanging the second component of the
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

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute<POD, 2>(*arr, f.on_device_);

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    // sum1
    POD sum1(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum1 += arrv_t(j, i) * arrv_t(j, i);
        }
    }

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 2>(*arr, f.on_device_);

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    // sum2
    POD sum2(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum2 += arrv_init_t(j, i) * arrv_t(j, i);
        }
    }

    adjoint_test(sum1, sum2, "test_rank1_strided_v2_adj_test");
}


void test_rank2(Fixture& f) {
    array::ArrayT<POD> arr(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv = array::make_host_view<POD, 3>(arr);


    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 3>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0,  0,    0,   0,  0, 0,   -2, 2,    -20, 20, -200, 200, -4, 4, -40,
                           40, -400, 400, -6, 6, -60, 60, -600, 600, 0,  0,    0,   0,  0, 0};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0,   0,  0,    0,   0,     0,    -8, 8, -80, 80, -800, 800, -10, 10, -100, 100, -1000, 1000,
                           -12, 12, -120, 120, -1200, 1200, 0,  0, 0,   0,  0,    0,   0,   0,  0,    0,   0,     0};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0,     0,    0,     0,    0,   0,  0,    0,   0,     0,    0,   0,  -14,  14,
                           -140,  140,  -1400, 1400, -16, 16, -160, 160, -1600, 1600, -18, 18, -180, 180,
                           -1800, 1800, 0,     0,    0,   0,  0,    0,   0,     0,    0,   0,  0,    0};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank2_adj_test(Fixture& f) {
    array::ArrayT<POD> arr_init(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_init = array::make_host_view<POD, 3>(arr_init);
    array::ArrayT<POD> arr(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv = array::make_host_view<POD, 3>(arr);

    for (size_t p = 0; p < static_cast<size_t>(f.N); ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv_init(p, i, static_cast<size_t>(0)) =
                (static_cast<size_t>(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv_init(p, i, static_cast<size_t>(1)) =
                (static_cast<size_t>(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
            arrv(p, i, static_cast<size_t>(0)) = arrv_init(p, i, static_cast<size_t>(0));
            arrv(p, i, static_cast<size_t>(1)) = arrv_init(p, i, static_cast<size_t>(1));
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute<POD, 3>(arr, f.on_device_);

    // sum1
    POD sum1(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 2; ++k) {
                sum1 += arrv(p, i, k) * arrv(p, i, k);
            }
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 3>(arr, f.on_device_);

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    // sum2
    POD sum2(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 2; ++k) {
                sum2 += arrv_init(p, i, k) * arrv(p, i, k);
            }
        }
    }

    adjoint_test(sum1, sum2, "test_rank2_adj_test");
}

void test_rank2_l1(Fixture& f) {
    array::ArrayT<POD> arr_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_t = array::make_host_view<POD, 3>(arr_t);
    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
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

    f.halo_exchange_std->execute_adjoint<POD, 3>(*arr, false);


    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0,  0, -90, 90, -900, 900,   // halo
                           -2, 2, -10, 10, -100, 100,   // core
                           -4, 4, -20, 20, -200, 200,   // core
                           -6, 6, -30, 30, -300, 300,   // core
                           0,  0, -40, 40, -400, 400};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0,   0,  -30, 30, -300, 300,   // halo
                           -8,  8,  -40, 40, -400, 400,   // core
                           -10, 10, -50, 50, -500, 500,   // core
                           -12, 12, -60, 60, -600, 600,   // core
                           0,   0,  -70, 70, -700, 700,   // halo
                           0,   0,  -80, 80, -800, 800};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0,   0,  -50, 50, -500, 500,   // halo
                           0,   0,  -60, 60, -600, 600,   // halo
                           -14, 14, -70, 70, -700, 700,   // core
                           -16, 16, -80, 80, -800, 800,   // core
                           -18, 18, -90, 90, -900, 900,   // core
                           0,   0,  -10, 10, -100, 100,   // halo
                           0,   0,  -20, 20, -200, 200};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
    }
}

void test_rank2_l1_adj_test(Fixture& f) {
    array::ArrayT<POD> arr_init_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_init_t = array::make_host_view<POD, 3>(arr_init_t);
    array::ArrayT<POD> arr_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_t = array::make_host_view<POD, 3>(arr_t);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 3; ++i) {
            arrv_init_t(p, i, static_cast<std::size_t>(0)) =
                (static_cast<std::size_t>(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv_init_t(p, i, static_cast<std::size_t>(1)) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
            arrv_t(p, i, static_cast<std::size_t>(0)) = arrv_init_t(p, i, static_cast<std::size_t>(0));
            arrv_t(p, i, static_cast<std::size_t>(1)) = arrv_init_t(p, i, static_cast<std::size_t>(1));
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

    f.halo_exchange_std->execute<POD, 3>(*arr, false);

    // sum1
    POD sum1(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 2; ++k) {
                sum1 += arrv_t(p, i, k) * arrv_t(p, i, k);
            }
        }
    }

    f.halo_exchange_std->execute_adjoint<POD, 3>(*arr, false);

    // sum2
    POD sum2(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 2; ++k) {
                sum2 += arrv_init_t(p, i, k) * arrv_t(p, i, k);
            }
        }
    }

    adjoint_test(sum1, sum2, "test_rank2_l1_adj_test");
}

void test_rank2_l2_v2(Fixture& f) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    // Test rank 2 halo-exchange
    array::ArrayT<POD> arr_t(f.N, 3, 2);
    array::ArrayView<POD, 3> arrv_t = array::make_host_view<POD, 3>(arr_t);
    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
    }

    /*
    for ( int p = 0; p < f.N; ++p ) {
        for ( size_t i = 0; i < 3; ++i ) {
            arrv_t( (size_t)p, i, (size_t)0 ) =
                ( size_t( f.part[p] ) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow( 10, i ) );
            arrv_t( (size_t)p, i, (size_t)1 ) =
                ( size_t( f.part[p] ) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow( 10, i ) );
        }
    }
   */

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        &arrv_t(0, 1, 1), array::ArraySpec {
            array::make_shape(f.N, 1, 1),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
                array::make_strides(192, 32, 1)
#else
                     array::make_strides(6, 2, 1)
#endif
        }));

    f.halo_exchange_std->execute_adjoint<POD, 3>(*arr, f.on_device_);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {-9, 9, -90, 0,  -900, 900,   // halo
                           -1, 1, -10, 20, -100, 100,   // core
                           -2, 2, -20, 40, -200, 200,   // core
                           -3, 3, -30, 60, -300, 300,   // core
                           -4, 4, -40, 0,  -400, 400};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {-3,
                           3,
                           -30,
                           0,
                           -300,
                           300,  // halo
                           -4,
                           4,
                           -40,
                           80,
                           -400,
                           400,  // core
                           -5,
                           5,
                           -50,
                           100,
                           -500,
                           500,  // core
                           -6,
                           6,
                           -60,
                           120,
                           -600,
                           600,  // core
                           -7,
                           7,
                           -70,
                           0,
                           -700,
                           700  // halo
                               - 8,
                           8,
                           -80,
                           0,
                           -800,
                           800};  // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {-5, 5, -50, 0,   -500, 500,   // halo
                           -6, 6, -60, 0,   -600, 600,   // halo
                           -7, 7, -70, 140, -700, 700,   // core
                           -8, 8, -80, 160, -800, 800,   // core
                           -9, 9, -90, 180, -900, 900,   // core
                           -1, 1, -10, 0,   -100, 100,   // halo
                           -2, 2, -20, 0,   -200, 200};  // halo
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

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv_t(j, i, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
    }
    /*
    for ( int p = 0; p < f.N; ++p ) {
        for ( size_t i = 0; i < 3; ++i ) {
            arrv_t( (size_t)p, i, (size_t)0 ) =
                ( size_t( f.part[p] ) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow( 10, i ) );
            arrv_t( (size_t)p, i, (size_t)1 ) =
                ( size_t( f.part[p] ) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow( 10, i ) );
        }
    }
*/

    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(
        &arrv_t(0, 0, 1), array::ArraySpec {
            array::make_shape(f.N, 3, 1),
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
                array::make_strides(192, 32, 2)
#else
                     array::make_strides(6, 2, 2)
#endif
        }));


    f.halo_exchange_std->execute_adjoint<POD, 3>(*arr, f.on_device_);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {-9, 0, -90, 0,  -900, 0,    // halo
                           -1, 2, -10, 20, -100, 200,  // core
                           -2, 4, -20, 40, -200, 400,  // core
                           -3, 6, -30, 60, -300, 600,  // core
                           -4, 0, -40, 0,  -400, 0};   // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {-3, 0,  -30, 0,   -300, 0,     // halo
                           -4, 8,  -40, 80,  -400, 800,   // core
                           -5, 10, -50, 100, -500, 1000,  // core
                           -6, 12, -60, 120, -600, 1200,  // core
                           -7, 0,  -70, 0,   -700, 0,     // halo
                           -8, 0,  -80, 0,   -800, 0};    // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {-5, 0,  -50, 0,   -500, 0,     // halo
                           -6, 0,  -60, 0,   -600, 0,     // halo
                           -7, 14, -70, 140, -700, 1400,  // core
                           -8, 16, -80, 160, -800, 1600,  // core
                           -9, 18, -90, 180, -900, 1800,  // core
                           -1, 0,  -10, 0,   -100, 0,     // halo
                           -2, 0,  -20, 0,   -200, 0};    // halo
            validate<POD, 3>::apply(arrv_t, arr_c);
            break;
        }
    }
#endif
}

void test_rank0_wrap(Fixture& f) {
    std::vector<POD> existing = f.gidx;
    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(existing.data(), array::make_shape(f.N)));
    array::ArrayView<POD, 1> arrv = array::make_view<POD, 1>(*arr);
    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 2, 4, 6, 0};
            for (int j = 0; j < f.N; ++j) {
                arrv(j) = arr_c[j];
            }
            break;
        }
        case 1: {
            POD arr_c[] = {0, 8, 10, 12, 0, 0};
            for (int j = 0; j < f.N; ++j) {
                arrv(j) = arr_c[j];
            }
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 14, 16, 18, 0, 0};
            for (int j = 0; j < f.N; ++j) {
                arrv(j) = arr_c[j];
            }
            break;
        }
    }

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 1>(*arr, f.on_device_);

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 2, 4, 6, 0};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 8, 10, 12, 0, 0};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 14, 16, 18, 0, 0};
            validate<POD, 1>::apply(arrv, arr_c);
            break;
        }
    }
    arr->deallocateDevice();
}

void test_rank0_wrap_adj_test(Fixture& f) {
    std::unique_ptr<array::Array> arr(array::Array::wrap<POD>(f.gidx.data(), array::make_shape(f.N)));
    array::ArrayView<POD, 1> arrv = array::make_view<POD, 1>(*arr);

    // note we have to be VERY careful here
    // we can't do
    //std::unique_ptr<array::Array> arr_init(
    //     array::Array::wrap<POD>( f.gidx.data(), array::make_shape( f.N ) ) );
    //array::ArrayView<POD, 1> arrv_init =
    //     array::make_view<POD, 1>( *arr_init );
    //
    // as both pointers will end up sharing the same memory !!
    //
    // but instead we need to do;

    array::ArrayT<POD> arr_init(f.N);
    array::ArrayView<POD, 1> arrv_init = array::make_host_view<POD, 1>(arr_init);

    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        arrv_init(j) = arrv(j);
    }

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute<POD, 1>(*arr, f.on_device_);

    // sum1
    POD sum1(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        sum1 += arrv(j) * arrv(j);
    }

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    f.halo_exchange_std->execute_adjoint<POD, 1>(*arr, f.on_device_);

    if (f.on_device_) {
        arr->syncHostDevice();
    }

    // sum2
    POD sum2(0);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        sum2 += arrv_init(j) * arrv(j);
    }

    adjoint_test(sum1, sum2, "test_rank0_wrap_adj_test");
}

void test_rank1_paralleldim1(Fixture& f) {
    array::ArrayT<POD> arr(2, f.N);
    array::ArrayView<POD, 2> arrv = array::make_view<POD, 2>(arr);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (int j = 0; j < f.N; ++j) {
                arrv(0, j) = arr_c[j] * 10;
                arrv(1, j) = arr_c[j] * 100;
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (int j = 0; j < f.N; ++j) {
                arrv(0, j) = arr_c[j] * 10;
                arrv(1, j) = arr_c[j] * 100;
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (int j = 0; j < f.N; ++j) {
                arrv(0, j) = arr_c[j] * 10;
                arrv(1, j) = arr_c[j] * 100;
            }
            break;
        }
    }

    f.halo_exchange_std->execute_adjoint<POD, 2, array::LastDim>(arr, false);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 20, 40, 60, 0, 0, 200, 400, 600, 0};  // 90,900, 10,100, 20,200, 30,300, 40,400 };
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 80, 100, 120, 0, 0, 0, 800, 1000, 1200, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 140, 160, 180, 0, 0, 0, 0, 1400, 1600, 1800, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank1_paralleldim1_adj_test(Fixture& f) {
    array::ArrayT<POD> arr_init(2, f.N);
    array::ArrayView<POD, 2> arrv_init = array::make_view<POD, 2>(arr_init);
    array::ArrayT<POD> arr(2, f.N);
    array::ArrayView<POD, 2> arrv = array::make_view<POD, 2>(arr);
    for (std::size_t j = 0; j < static_cast<std::size_t>(f.N); ++j) {
        arrv_init(0ul, j) = (static_cast<std::size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 10);
        arrv_init(1ul, j) = (static_cast<std::size_t>(f.part[j]) != mpi::comm().rank() ? 0 : f.gidx[j] * 100);
        arrv(0ul, j)      = arrv_init(0ul, j);
        arrv(1ul, j)      = arrv_init(1ul, j);
    }

    f.halo_exchange_std->execute<POD, 2, array::LastDim>(arr, false);

    // sum1
    POD sum1(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum1 += arrv(i, p) * arrv(i, p);
        }
    }

    f.halo_exchange_std->execute_adjoint<POD, 2, array::LastDim>(arr, false);

    // sum2
    POD sum2(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 2; ++i) {
            sum2 += arrv_init(i, p) * arrv(i, p);
        }
    }

    adjoint_test(sum1, sum2, "test_rank1_paralleldim1_adj_test");
}

void test_rank2_paralleldim2(Fixture& f) {
    array::ArrayT<POD> arr(3, f.N, 2);
    array::ArrayView<POD, 3> arrv = array::make_view<POD, 3>(arr);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv(i, j, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv(i, j, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (std::size_t j = 0; j < static_cast<size_t>(f.N); ++j) {
                for (std::size_t i = 0; i < 3; ++i) {
                    for (std::size_t k = 0; k < 2; ++k) {
                        arrv(i, j, k) = (2 * static_cast<int>(k) - 1) * arr_c[j] * std::pow(10, i);
                    }
                }
            }
            break;
        }
    }

    f.halo_exchange_std->execute_adjoint<POD, 3, array::Dim<1>>(arr, false);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0,  0,   -2, 2, -4, 4, -6, 6,    0,   0,    0,   0,    -20, 20, -40,
                           40, -60, 60, 0, 0,  0, 0,  -200, 200, -400, 400, -600, 600, 0,  0};

            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0,    0,   -8, 8, -10, 10, -12, 12, 0,    0,   0,     0,    0,     0,    -80, 80, -100, 100,
                           -120, 120, 0,  0, 0,   0,  0,   0,  -800, 800, -1000, 1000, -1200, 1200, 0,   0,  0,    0};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 0, 0, -14,   14,   -16,   16,   -18,   18,   0, 0, 0, 0,
                           0, 0, 0, 0, -140,  140,  -160,  160,  -180,  180,  0, 0, 0, 0,
                           0, 0, 0, 0, -1400, 1400, -1600, 1600, -1800, 1800, 0, 0, 0, 0};
            validate<POD, 3>::apply(arrv, arr_c);
            break;
        }
    }
}

void test_rank2_paralleldim2_adj_test(Fixture& f) {
    array::ArrayT<POD> arr_init(3, f.N, 2);
    array::ArrayView<POD, 3> arrv_init = array::make_view<POD, 3>(arr_init);
    array::ArrayT<POD> arr(3, f.N, 2);
    array::ArrayView<POD, 3> arrv = array::make_view<POD, 3>(arr);
    for (size_t p = 0; p < static_cast<size_t>(f.N); ++p) {
        for (size_t i = 0; i < 3; ++i) {
            arrv_init(i, p, static_cast<size_t>(0)) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : -f.gidx[p] * std::pow(10, i));
            arrv_init(i, p, static_cast<size_t>(1)) =
                (size_t(f.part[p]) != mpi::comm().rank() ? 0 : f.gidx[p] * std::pow(10, i));
            arrv(i, p, static_cast<size_t>(0)) = arrv_init(i, p, static_cast<size_t>(0));
            arrv(i, p, static_cast<size_t>(1)) = arrv_init(i, p, static_cast<size_t>(1));
        }
    }

    f.halo_exchange_std->execute<POD, 3, array::Dim<1>>(arr, false);


    // sum1
    POD sum1(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 2; ++k) {
                sum1 += arrv(i, p, k) * arrv(i, p, k);
            }
        }
    }

    f.halo_exchange_std->execute_adjoint<POD, 3, array::Dim<1>>(arr, false);

    // sum2
    POD sum2(0);
    for (std::size_t p = 0; p < static_cast<std::size_t>(f.N); ++p) {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t k = 0; k < 2; ++k) {
                sum2 += arrv_init(i, p, k) * arrv(i, p, k);
            }
        }
    }

    adjoint_test(sum1, sum2, "test_rank2_paralleldim2_adj_test");
}


void test_rank1_cinterface(Fixture& f) {
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST

    std::cout << "entering test_rank1_cinterface" << std::endl;

    array::ArrayT<POD> arr(f.N, 2);

    array::ArrayView<POD, 2> arrv = array::make_host_view<POD, 2>(arr);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {9, 1, 2, 3, 4};
            for (int j = 0; j < f.N; ++j) {
                arrv(j, 0) = arr_c[j] * 10;
                arrv(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 1: {
            POD arr_c[] = {3, 4, 5, 6, 7, 8};
            for (int j = 0; j < f.N; ++j) {
                arrv(j, 0) = arr_c[j] * 10;
                arrv(j, 1) = arr_c[j] * 100;
            }
            break;
        }
        case 2: {
            POD arr_c[] = {5, 6, 7, 8, 9, 1, 2};
            for (int j = 0; j < f.N; ++j) {
                arrv(j, 0) = arr_c[j] * 10;
                arrv(j, 1) = arr_c[j] * 100;
            }
            break;
        }
    }

    if (f.on_device_) {
        arr.syncHostDevice();
    }

    int shapes[2]  = {(int)arrv.shape(0), (int)arrv.shape(1)};
    int strides[2] = {(int)arrv.stride(0), (int)arrv.stride(1)};

    atlas__HaloExchange__execute_adjoint_strided_double(f.halo_exchange_std.get(), arrv.data(), &(strides[1]),
                                                        &(shapes[1]), 1);

    switch (mpi::comm().rank()) {
        case 0: {
            POD arr_c[] = {0, 0, 20, 200, 40, 400, 60, 600, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 1: {
            POD arr_c[] = {0, 0, 80, 800, 100, 1000, 120, 1200, 0, 0, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
        case 2: {
            POD arr_c[] = {0, 0, 0, 0, 140, 1400, 160, 1600, 180, 1800, 0, 0, 0, 0};
            validate<POD, 2>::apply(arrv, arr_c);
            break;
        }
    }
#endif
}

CASE("test_haloexchange_adjoint") {
    Fixture f(false);

    SECTION("test_rank0_arrview") { test_rank0_arrview(f); }

    SECTION("test_rank0_arrview_adj_test") { test_rank0_arrview_adj_test(f); }

    SECTION("test_rank1") { test_rank1(f); }

    SECTION("test_rank1_adj_test") { test_rank1_adj_test(f); }

    SECTION("test_rank1_strided_v1") { test_rank1_strided_v1(f); }

    SECTION("test_rank1_strided_v1_adj_test") { test_rank1_strided_v1_adj_test(f); }

    SECTION("test_rank1_strided_v2") { test_rank1_strided_v2(f); }

    SECTION("test_rank1_strided_v2_adj_test") { test_rank1_strided_v2_adj_test(f); }

    SECTION("test_rank2") { test_rank2(f); }

    SECTION("test_rank2_adj_test") { test_rank2_adj_test(f); }

    SECTION("test_rank2_l1") { test_rank2_l1(f); }

    SECTION("test_rank2_l1_adj_test") { test_rank2_l1_adj_test(f); }

    /* NOT TESTED AS USES ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    SECTION( "test_rank2_l2_v2" ) { test_rank2_l2_v2( f ); }

    SECTION( "test_rank2_v2" ) { test_rank2_v2( f ); }
*/

    SECTION("test_rank0_wrap") { test_rank0_wrap(f); }

    SECTION("test_rank0_wrap_adj_test") { test_rank0_wrap_adj_test(f); }

    SECTION("test_rank1_paralleldim_1") { test_rank1_paralleldim1(f); }

    SECTION("test_rank1_paralleldim_1_adj_test") { test_rank1_paralleldim1_adj_test(f); }

    SECTION("test_rank2_paralleldim_2") { test_rank2_paralleldim2(f); }

    SECTION("test_rank2_paralleldim_2_adj_test") { test_rank2_paralleldim2_adj_test(f); }

    /* NOT TESTED AS USES ATLAS_GRIDTOOLS_STORAGE_BACKEND_HOST
    SECTION( "test_rank1_cinterface" ) { test_rank1_cinterface( f ); }
*/

    /* NOT IMPLEMENTED
#if ATLAS_GRIDTOOLS_STORAGE_BACKEND_CUDA
    f.on_device_ = true;

    SECTION( "test_rank0_arrview" ) { test_rank0_arrview( f ); }

    SECTION( "test_rank1" ) { test_rank1( f ); }

    SECTION( "test_rank2" ) { test_rank2( f ); }
    SECTION( "test_rank0_wrap" ) { test_rank0_wrap( f ); }

#endif
*/
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
