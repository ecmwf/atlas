/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/runtime/Exception.h"
#include "atlas/array.h"
#include "hic/hic.h"

namespace atlas {
namespace detail {

#define USE_MDSPAN 0
template <class Blocked, class Nonblocked>
__global__ void kernel_copy_blocked_to_nonblocked_mdspan(const Blocked blocked, Nonblocked nonblocked) {
    auto np     = nonblocked.extent(0);
    auto nblks  = blocked.extent(0);
    auto nproma = blocked.extent(blocked.rank()-1);
    static_assert(nonblocked.rank() == blocked.rank()-1);
    if constexpr(blocked.rank()==4) {
        // ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(2));
        // ATLAS_ASSERT(nonblocked.extent(2) == blocked.extent(1));
        idx_t nlev = nonblocked.extent(1);
        idx_t nvar = nonblocked.extent(2);
        for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
            idx_t nrof;
            if (np - jpbegin < nproma) {
                nrof = np - jpbegin;
            }
            else {
                nrof = nproma;
            }
            // auto nrof = std::min(np - jpbegin, nproma);
            for (idx_t jvar = 0; jvar < nvar; ++jvar) {
                for (idx_t jlev = 0; jlev < nlev; ++jlev) {
                    for (idx_t jrof = 0; jrof < nrof; ++jrof) {
                        auto jp = jpbegin+jrof;
                        nonblocked(jp, jlev, jvar) = blocked(jblk, jvar, jlev, jrof);
                    }
                }
            }
        }
    }
    // else if constexpr (blocked.rank()==3) {
    //     ATLAS_ASSERT(nonblocked.extent(1) == blocked.extent(1));
    //     idx_t nlev = nonblocked.extent(1);
    //     for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
    //         auto nrof = std::min(np - jpbegin, nproma);
    //         for (idx_t jlev = 0; jlev < nlev; ++jlev) {
    //             for (idx_t jrof = 0; jrof < nrof; ++jrof) {
    //                 auto jp = jpbegin+jrof;
    //                 nonblocked(jp, jlev) = blocked(jblk, jlev, jrof);
    //             }
    //         }
    //     }
    // }
    // else if constexpr (blocked.rank()==2) {
    //     for (idx_t jblk = 0, jpbegin = 0; jblk < nblks; ++jblk, jpbegin+=nproma) {
    //         auto nrof = std::min(np - jpbegin, nproma);
    //         for (idx_t jrof = 0; jrof < nrof; ++jrof) {
    //             auto jp = jpbegin+jrof;
    //             nonblocked(jp) = blocked(jblk, jrof);
    //         }
    //     }
    // }
    else {
        // ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}


// template<int ParallelDim, typename DATA_TYPE, int RANK>
// __global__ void pack_kernel(const int sendcnt, const int* sendmap_ptr, const idx_t sendmap_size,
//          const array::ArrayView<DATA_TYPE, RANK> field, DATA_TYPE* send_buffer,
//          const idx_t send_buffer_size, const typename std::enable_if<RANK==1, int>::type = 0) {
//     const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);

//     const idx_t node_cnt = blockIdx.x*blockDim.x + threadIdx.x;

//     if(node_cnt >= sendcnt) return;

//     send_buffer[node_cnt] = field(sendmap[node_cnt]);
// }

// template<int ParallelDim, typename DATA_TYPE, int RANK>
// __global__ void pack_kernel(const int sendcnt,  const int* sendmap_ptr, const idx_t sendmap_size,
//          const array::ArrayView<DATA_TYPE, RANK> field, DATA_TYPE* send_buffer,
//          const idx_t send_buffer_size, const typename std::enable_if<RANK>=2, int>::type = 0) {
//     const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);

//     const idx_t node_cnt = blockIdx.x*blockDim.x + threadIdx.x;
//     const idx_t var1_idx = blockIdx.y*blockDim.y + threadIdx.y;

//     if(node_cnt >= sendcnt || var1_idx >= field.shape(1) ) return;

//     idx_t buff_idx = get_buffer_index<ParallelDim, RANK>::apply(field, node_cnt, var1_idx);
//     const idx_t node_idx = sendmap[node_cnt];

//     pack_index<0, (RANK>=3) ? (RANK-2) : 0, 2>::apply(buff_idx, node_idx, field, send_buffer, node_idx,var1_idx);
// }

// template<typename Blocked, typename NonBlocked>
// struct Operator {
//     void apply( Blocked blocked, NonBlocked nonblocked) {

//     const unsigned int block_size_x = 32;
//     const unsigned int block_size_y = (RANK==1) ? 1 : 4;

//     const unsigned int nblocks_x = (sendcnt+block_size_x-1)/block_size_x;
//     const unsigned int nblocks_y = get_n_hic_blocks<ParallelDim, RANK>::apply(array, block_size_y);

//     dim3 threads(block_size_x, block_size_y);
//     dim3 blocks(nblocks_x, nblocks_y);
//     HIC_CALL(hicDeviceSynchronize());
//     hicError_t err = hicGetLastError();
//     if (err != hicSuccess) {
//         std::string msg = std::string("Error synchronizing device")+ hicGetErrorString(err);
//         throw_Exception(msg);
//     }

//   pack_kernel<ParallelDim, DATA_TYPE, RANK><<<blocks,threads>>>(sendcnt, sendmap.data(), sendmap.size(), array, send_buffer, send_buffer_size);
//     err = hicGetLastError();
//     if (err != hicSuccess)
//         throw_Exception("Error launching GPU packing kernel");

//     HIC_CALL(hicDeviceSynchronize());
//     err = hicGetLastError();
//     if (err != hicSuccess) {
//         std::string msg = std::string("Error synchronizing device")+ hicGetErrorString(err);
//         throw_Exception(msg);
//     }
// };


template <class ValueType>
void copy_blocked_to_nonblocked_T(const array::Array& blocked, array::Array& nonblocked, bool on_device) {
    ATLAS_ASSERT(nonblocked.rank() == blocked.rank()-1);
    if (blocked.rank()==4) {
        auto blocked_v    = on_device ? array::make_device_view<ValueType, 4>(blocked)    : array::make_host_view<ValueType, 4>(blocked);
        auto nonblocked_v = on_device ? array::make_device_view<ValueType, 3>(nonblocked) : array::make_host_view<ValueType, 3>(nonblocked);
        if(USE_MDSPAN) {
            // copy_blocked_to_nonblocked_mdspan(blocked_v.as_mdspan(), nonblocked_v.as_mdspan(), on_device);
        }
        else {
            kernel_copy_blocked_to_nonblocked_mdspan<<<1,1>>>(blocked_v, nonblocked_v);
        }
    }
    // else if (blocked.rank()==3) {
    //     auto blocked_v    = on_device ? array::make_device_view<ValueType, 3>(blocked)    : array::make_host_view<ValueType, 3>(blocked);
    //     auto nonblocked_v = on_device ? array::make_device_view<ValueType, 2>(nonblocked) : array::make_host_view<ValueType, 2>(nonblocked);
    //     if(USE_MDSPAN) {
    //         copy_blocked_to_nonblocked_mdspan(blocked_v.as_mdspan(), nonblocked_v.as_mdspan(), on_device);
    //     }
    //     else {
    //         copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
    //     }
    // }
    // else if (blocked.rank()==2) {
    //     auto blocked_v    = on_device ? array::make_device_view<ValueType, 2>(blocked)    : array::make_host_view<ValueType, 2>(blocked);
    //     auto nonblocked_v = on_device ? array::make_device_view<ValueType, 1>(nonblocked) : array::make_host_view<ValueType, 1>(nonblocked);
    //     if(USE_MDSPAN) {
    //         copy_blocked_to_nonblocked_mdspan(blocked_v.as_mdspan(), nonblocked_v.as_mdspan(), on_device);
    //     }
    //     else {
    //         copy_blocked_to_nonblocked_mdspan(blocked_v, nonblocked_v, on_device);
    //     }
    // }
    else {
        ATLAS_THROW_EXCEPTION("transposition not implemented");
    }
}

void copy_blocked_to_nonblocked(const array::Array& blocked, array::Array& nonblocked, bool on_device) {
    ATLAS_ASSERT(blocked.datatype() == nonblocked.datatype());
    switch (blocked.datatype().kind()) {
        // case array::DataType::kind<int>()    : return copy_blocked_to_nonblocked_T<int>(blocked, nonblocked, on_device);
        // case array::DataType::kind<long>()   : return copy_blocked_to_nonblocked_T<long>(blocked, nonblocked, on_device);
        // case array::DataType::kind<float>()  : return copy_blocked_to_nonblocked_T<float>(blocked, nonblocked, on_device);
        case array::DataType::kind<double>() : return copy_blocked_to_nonblocked_T<double>(blocked, nonblocked, on_device);
        default: throw_Exception("datatype not supported", Here());
    }
}

}
} //namespace atlas::parallel::detail
