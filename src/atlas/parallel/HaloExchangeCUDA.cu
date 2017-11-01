/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "HaloExchangeCUDA.h"

namespace atlas {
namespace parallel {

template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt, const array::SVector<int> sendmap,
         const array::ArrayView<DATA_TYPE, RANK, false> field, array::SVector<DATA_TYPE> send_buffer,
                            typename std::enable_if<RANK==2, int>::type* = 0) {
    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= sendcnt || i >= field.data_view().template length<1>() ) return;

    const size_t buff_idx = field.data_view().template length<1>() * p + i;

    send_buffer[buff_idx] = field(sendmap[p], i);
}

template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt, const array::SVector<int> sendmap,
         const array::ArrayView<DATA_TYPE, RANK, false> field, array::SVector<DATA_TYPE> send_buffer,
                            typename std::enable_if<RANK!=2, int>::type* = 0) {
}

template<typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int sendcnt, const array::SVector<int> recvmap,
         const array::SVector<DATA_TYPE> recv_buffer, array::ArrayView<DATA_TYPE, RANK> field,
                            typename std::enable_if<RANK==2, int>::type* = 0) {

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= sendcnt || i >= field.data_view().template length<1>() ) return;

    const size_t buff_idx = field.data_view().template length<1>() * p + i;

    field(recvmap[p], i) = recv_buffer[buff_idx];
}

template<typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int sendcnt, const array::SVector<int> recvmap,
         const array::SVector<DATA_TYPE> recv_buffer, array::ArrayView<DATA_TYPE, RANK> field,
                            typename std::enable_if<RANK!=2, int>::type* = 0) {
}

template<typename DATA_TYPE>
void halo_packer_cuda<DATA_TYPE, 1>::pack( const int sendcnt, array::SVector<int> const & sendmap,
                   const array::ArrayView<DATA_TYPE, 1>& field, array::SVector<DATA_TYPE>& send_buffer )
{
}

template<typename DATA_TYPE>
void halo_packer_cuda<DATA_TYPE, 1>::unpack(const int sendcnt, array::SVector<int> const & recvmap,
                   const array::SVector<DATA_TYPE> &recv_buffer ,
                   array::ArrayView<DATA_TYPE, 1> &field)
{
}

template<typename DATA_TYPE, int RANK>
void halo_packer_cuda<DATA_TYPE, RANK>::pack( const int sendcnt, array::SVector<int> const & sendmap,
                   const array::ArrayView<DATA_TYPE, RANK>& field, array::SVector<DATA_TYPE>& send_buffer )
{
  const unsigned int block_size_x = 32;
  const unsigned int block_size_y = 4;
  dim3 threads(block_size_x, block_size_y);
  dim3 blocks((sendcnt+block_size_x-1)/block_size_x, (field.data_view().template length<1>()+block_size_y-1)/block_size_y);

  pack_kernel<DATA_TYPE, RANK><<<blocks,threads>>>(sendcnt, sendmap, field, send_buffer);
}

template<typename DATA_TYPE, int RANK>
void halo_packer_cuda<DATA_TYPE, RANK>::unpack(const int sendcnt, array::SVector<int> const & recvmap,
                   const array::SVector<DATA_TYPE> &recv_buffer ,
                   array::ArrayView<DATA_TYPE, RANK> &field)
{
  const unsigned int block_size_x = 32;
  const unsigned int block_size_y = 4;
  dim3 threads(block_size_x, block_size_y);
  dim3 blocks((sendcnt+block_size_x-1)/block_size_x, (field.data_view().template length<1>()+block_size_y-1)/block_size_y);

  unpack_kernel<<<blocks,threads>>>(sendcnt, recvmap, recv_buffer, field);
}

#define EXPLICIT_TEMPLATE_INSTANTIATION(RANK) \
template class halo_packer_cuda<int,RANK>; \
template class halo_packer_cuda<long,RANK>; \
template class halo_packer_cuda<long unsigned,RANK>; \
template class halo_packer_cuda<float,RANK>; \
template class halo_packer_cuda<double,RANK>; \

  EXPLICIT_TEMPLATE_INSTANTIATION(1)
  EXPLICIT_TEMPLATE_INSTANTIATION(2)
  EXPLICIT_TEMPLATE_INSTANTIATION(3)
  EXPLICIT_TEMPLATE_INSTANTIATION(4)

} //namespace array
} //namespace atlas
