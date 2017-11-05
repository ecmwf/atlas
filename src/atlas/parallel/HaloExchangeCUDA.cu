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
#include "atlas/array/SVector.h"

namespace atlas {
namespace parallel {

template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt,  const int* sendmap_ptr, const size_t sendmap_size,
         const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field, DATA_TYPE* send_buffer, const typename std::enable_if<RANK==2, int>::type = 0) {
    const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);
    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= sendcnt || i >= field.data_view().template length<1>() ) return;

    const size_t buff_idx = field.data_view().template length<1>() * p + i;

    send_buffer[buff_idx] = field(sendmap[p], i);

}

template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt,  const int* sendmap_ptr, const size_t sendmap_size,
         const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field, DATA_TYPE* send_buffer, const typename std::enable_if<RANK==3, int>::type = 0) {
    const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= sendcnt || i >= field.data_view().template length<1>() ) return;

    size_t buff_idx = field.data_view().template length<2>() * field.data_view().template length<1>() * p + field.data_view().template length<2>() * i;

    for(size_t varid=0; varid < field.data_view().template length<2>(); ++varid) {
        send_buffer[buff_idx++] = field(sendmap[p], i, varid);
    }

}

template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt, const int* sendmap_ptr, const size_t sendmap_size,
         const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field, DATA_TYPE* send_buffer,
                            const typename std::enable_if<RANK==1, int>::type = 0) {
    const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;

    if(p >= sendcnt) return;
    send_buffer[p] = field(sendmap[p]);
}

template<typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt, const int* sendmap, const size_t sendmap_size,
         const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field, DATA_TYPE* send_buffer,
                            const typename std::enable_if<(RANK>3), int>::type = 0) {
}

template<typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int recvcnt, const int* recvmap_ptr, const size_t recvmap_size,
         const DATA_TYPE* recv_buffer, array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field,
                            const typename std::enable_if<RANK==2, int>::type = 0) {

    const array::SVector<int> recvmap(const_cast<int*>(recvmap_ptr), recvmap_size);

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= recvcnt || i >= field.data_view().template length<1>() ) return;

    const size_t buff_idx = field.data_view().template length<1>() * p + i;

    field(recvmap[p], i) = recv_buffer[buff_idx];

}

template<typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int recvcnt, const int* recvmap_ptr, const size_t recvmap_size,
         const DATA_TYPE* recv_buffer, array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field,
                            const typename std::enable_if<RANK==3, int>::type = 0) {

    const array::SVector<int> recvmap(const_cast<int*>(recvmap_ptr), recvmap_size);

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;
    const size_t i = blockIdx.y*blockDim.y + threadIdx.y;

    if(p >= recvcnt || i >= field.data_view().template length<1>() ) return;

    size_t buff_idx = field.data_view().template length<2>() * field.data_view().template length<1>() * p + field.data_view().template length<2>() * i;

    for(size_t varid=0; varid < field.data_view().template length<2>(); ++varid) {
        field(recvmap[p], i, varid) = recv_buffer[buff_idx++];
    }
}

template<typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int recvcnt, const int* recvmap_ptr, const size_t recvmap_size,
         const DATA_TYPE* recv_buffer, array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field,
                            const typename std::enable_if<RANK==1, int>::type = 0) {

    const array::SVector<int> recvmap(const_cast<int*>(recvmap_ptr), recvmap_size);

    const size_t p = blockIdx.x*blockDim.x + threadIdx.x;

    if(p >= recvcnt) return;
    field(recvmap[p]) = recv_buffer[p];
}

template<typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int recvcnt, const int* recvmap_ptr, const size_t recvmap_size,
         const DATA_TYPE* recv_buffer, array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadWrite> field,
                            typename std::enable_if<(RANK>3), int>::type* = 0) {
}

template<int RANK>
struct get_n_cuda_blocks
{
  template<typename DATA_TYPE>
  static unsigned int apply(const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield, const unsigned int block_size_y) {
      return (hfield.data_view().template length<1>()+block_size_y-1)/block_size_y;
  }
};

template<>
struct get_n_cuda_blocks<1> {
    template<typename DATA_TYPE>
    static unsigned int apply(const array::ArrayView<DATA_TYPE, 1, array::Intent::ReadOnly>& hfield, const unsigned int block_size_y) {
        return 1;
    }
};

template<typename DATA_TYPE, int RANK>
void halo_packer_cuda<DATA_TYPE, RANK>::pack( const int sendcnt, array::SVector<int> const & sendmap,
                   const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield, const array::ArrayView<DATA_TYPE, RANK>& dfield,
                   array::SVector<DATA_TYPE>& send_buffer )
{
  const unsigned int block_size_x = 32;
  const unsigned int block_size_y = (RANK==1) ? 1 : 4;

  unsigned int nblocks_y = get_n_cuda_blocks<RANK>::apply(hfield, block_size_y);

  dim3 threads(block_size_x, block_size_y);
  dim3 blocks((sendcnt+block_size_x-1)/block_size_x, nblocks_y);
  cudaDeviceSynchronize();

  pack_kernel<DATA_TYPE, RANK><<<blocks,threads>>>(sendcnt, sendmap.data(), sendmap.size(), dfield, (send_buffer.data()));

  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) 
    throw eckit::Exception("Error launching GPU packing kernel");

}

template<typename DATA_TYPE, int RANK>
void halo_packer_cuda<DATA_TYPE, RANK>::unpack(const int recvcnt, array::SVector<int> const & recvmap,
                   const array::SVector<DATA_TYPE> &recv_buffer ,
                   const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly> &hfield, array::ArrayView<DATA_TYPE, RANK> &dfield)
{
  const unsigned int block_size_x = 32;
  const unsigned int block_size_y = (RANK==1) ? 1 : 4;

  unsigned int nblocks_y = get_n_cuda_blocks<RANK>::apply(hfield, block_size_y);

  dim3 threads(block_size_x, block_size_y);
  dim3 blocks((recvcnt+block_size_x-1)/block_size_x, nblocks_y);

  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error synchronizing device")+ cudaGetErrorString(err);
    throw eckit::Exception(msg);
  }

  unpack_kernel<<<blocks,threads>>>(recvcnt, recvmap.data(), recvmap.size(), recv_buffer.data(), dfield);

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error launching GPU packing kernel")+ cudaGetErrorString(err);
    throw eckit::Exception(msg);
  }

  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error synchronizing device")+ cudaGetErrorString(err);
    throw eckit::Exception(msg);
  }
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
