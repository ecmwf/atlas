/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/parallel/HaloExchangeCUDA.h"

namespace atlas {
namespace parallel {

#ifdef __CUDACC__
template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(int sendcnt, array::SVector<int> sendmap,
         const array::ArrayView<DATA_TYPE, RANK> field, array::SVector<DATA_TYPE> send_buffer) {

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= sendcnt || i >= field.data_view().template length<1>() ) return;

    const size_t buff_idx = field.data_view().template length<1>() * p + i;

    send_buffer[buff_idx] = field.data(p, i);
}

#endif

} // namespace parallel
} // namespace atlas
