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
#include "atlas/parallel/HaloExchangeImpl.h"
#include "atlas/array/SVector.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace parallel {

template<int ParallelDim, int RANK>
struct get_buffer_index{
    template<typename DATA_TYPE>
    ATLAS_HOST_DEVICE
    static idx_t apply(const array::ArrayView<DATA_TYPE, RANK>& field,
                           const idx_t node_cnt, const idx_t var1_idx) {
        return field.shape(RANK-1) * field.shape(RANK-2) * node_cnt +
                field.shape(RANK-1) * var1_idx;
    }
};

template<int ParallelDim>
struct get_buffer_index<ParallelDim, 2>{
    template<typename DATA_TYPE>
    ATLAS_HOST_DEVICE
    static idx_t apply(const array::ArrayView<DATA_TYPE, 2>& field,
                           const idx_t node_cnt, const idx_t var1_idx) {
        return field.shape(1) * node_cnt + var1_idx;
    }
};

template<int ParallelDim, typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt, const int* sendmap_ptr, const idx_t sendmap_size,
         const array::ArrayView<DATA_TYPE, RANK> field, DATA_TYPE* send_buffer,
         const idx_t send_buffer_size, const typename std::enable_if<RANK==1, int>::type = 0) {
    const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);

    const idx_t node_cnt = blockIdx.x*blockDim.x + threadIdx.x;

    if(node_cnt >= sendcnt) return;

    send_buffer[node_cnt] = field(sendmap[node_cnt]);
}

template<int ParallelDim, typename DATA_TYPE, int RANK>
__global__ void pack_kernel(const int sendcnt,  const int* sendmap_ptr, const idx_t sendmap_size,
         const array::ArrayView<DATA_TYPE, RANK> field, DATA_TYPE* send_buffer,
         const idx_t send_buffer_size, const typename std::enable_if<RANK>=2, int>::type = 0) {
    const array::SVector<int> sendmap(const_cast<int*>(sendmap_ptr), sendmap_size);

    const idx_t node_cnt = blockIdx.x*blockDim.x + threadIdx.x;
    const idx_t var1_idx = blockIdx.y*blockDim.y + threadIdx.y;

    if(node_cnt >= sendcnt || var1_idx >= field.shape(1) ) return;

    idx_t buff_idx = get_buffer_index<ParallelDim, RANK>::apply(field, node_cnt, var1_idx);
    const idx_t node_idx = sendmap[node_cnt];

    halo_packer_impl<0, (RANK>=3) ? (RANK-2) : 0, 2>::apply(buff_idx, node_idx, field, send_buffer, node_idx,var1_idx);
}

template<int ParallelDim, typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int recvcnt, const int* recvmap_ptr, const idx_t recvmap_size,
         const DATA_TYPE* recv_buffer, const idx_t recv_buffer_size, array::ArrayView<DATA_TYPE, RANK> field, const typename std::enable_if<RANK==1, int>::type = 0) {

    const array::SVector<int> recvmap(const_cast<int*>(recvmap_ptr), recvmap_size);

    idx_t node_cnt = blockIdx.x*blockDim.x + threadIdx.x;

    if(node_cnt >= recvcnt) return;

    const idx_t node_idx = recvmap[node_cnt];

    field(node_idx) = recv_buffer[node_cnt];
}

template<int ParallelDim, typename DATA_TYPE, int RANK>
__global__ void unpack_kernel(const int recvcnt, const int* recvmap_ptr, const idx_t recvmap_size,
         const DATA_TYPE* recv_buffer, const idx_t recv_buffer_size, array::ArrayView<DATA_TYPE, RANK> field, const typename std::enable_if<RANK>=2, int>::type = 0) {

    const array::SVector<int> recvmap(const_cast<int*>(recvmap_ptr), recvmap_size);

    const idx_t node_cnt = blockIdx.x*blockDim.x + threadIdx.x;
    const idx_t var1_idx = blockIdx.y*blockDim.y + threadIdx.y;

    if(node_cnt >= recvcnt || var1_idx >= field.shape(1) ) return;

    const idx_t node_idx = recvmap[node_cnt];

    idx_t buff_idx = get_buffer_index<ParallelDim, RANK>::apply(field, node_cnt, var1_idx);

    halo_unpacker_impl<0, (RANK>=3) ? (RANK-2) : 0, 2>::apply(buff_idx, node_idx, recv_buffer, field,node_idx,var1_idx);
}

template<int ParallelDim, int RANK, int DimCnt>
struct get_first_non_parallel_dim
{
    static_assert((ParallelDim <= RANK), "Error: parallelDim larger than RANK");
    constexpr static int apply() {
        return (DimCnt == ParallelDim) ? get_first_non_parallel_dim<ParallelDim, RANK, DimCnt+1>::apply() : DimCnt;
    }
};

template<int ParallelDim, int RANK>
struct get_first_non_parallel_dim<ParallelDim, RANK, RANK>
{
    static_assert((ParallelDim <= RANK), "Error: parallelDim larger than RANK");
    constexpr static int apply() {
        return -1;
    }
};

template<int ParallelDim, int RANK>
struct get_n_cuda_blocks
{
  template<typename DATA_TYPE>
  static unsigned int apply(const array::ArrayView<DATA_TYPE, RANK>& hfield, const unsigned int block_size_y) {
      constexpr auto dim = get_first_non_parallel_dim<ParallelDim, RANK, 0>::apply();
      return (hfield.shape(dim)+block_size_y-1)/block_size_y;
  }
};

template<>
struct get_n_cuda_blocks<0, 1> {
    template<typename DATA_TYPE>
    static unsigned int apply(const array::ArrayView<DATA_TYPE, 1>& hfield, const unsigned int block_size_y) {
        return 1;
    }
};

template<int ParallelDim, typename DATA_TYPE, int RANK>
void halo_packer_cuda<ParallelDim, DATA_TYPE, RANK>::pack( const int sendcnt, array::SVector<int> const & sendmap,
                   const array::ArrayView<DATA_TYPE, RANK>& hfield, const array::ArrayView<DATA_TYPE, RANK>& dfield,
                   DATA_TYPE* send_buffer, int send_buffer_size )
{
  const unsigned int block_size_x = 32;
  const unsigned int block_size_y = (RANK==1) ? 1 : 4;

  unsigned int nblocks_y = get_n_cuda_blocks<ParallelDim, RANK>::apply(hfield, block_size_y);

  dim3 threads(block_size_x, block_size_y);
  dim3 blocks((sendcnt+block_size_x-1)/block_size_x, nblocks_y);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error synchronizing device")+ cudaGetErrorString(err);
    throw_Exception(msg);
  }

  pack_kernel<ParallelDim, DATA_TYPE, RANK><<<blocks,threads>>>(sendcnt, sendmap.data(), sendmap.size(), dfield, send_buffer, send_buffer_size);
  err = cudaGetLastError();
  if (err != cudaSuccess)
    throw_Exception("Error launching GPU packing kernel");

  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error synchronizing device")+ cudaGetErrorString(err);
    throw_Exception(msg);
  }

}

template<int ParallelDim, typename DATA_TYPE, int RANK>
void halo_packer_cuda<ParallelDim, DATA_TYPE, RANK>::unpack(const int recvcnt, array::SVector<int> const & recvmap,
                   const DATA_TYPE* recv_buffer , int recv_buffer_size,
                   array::ArrayView<DATA_TYPE, RANK> &hfield, array::ArrayView<DATA_TYPE, RANK> &dfield)
{
  const unsigned int block_size_x = 32;
  const unsigned int block_size_y = (RANK==1) ? 1 : 4;

  unsigned int nblocks_y = get_n_cuda_blocks<ParallelDim, RANK>::apply(hfield, block_size_y);

  dim3 threads(block_size_x, block_size_y);
  dim3 blocks((recvcnt+block_size_x-1)/block_size_x, nblocks_y);

  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error synchronizing device")+ cudaGetErrorString(err);
    throw_Exception(msg);
  }

  unpack_kernel<ParallelDim, DATA_TYPE, RANK><<<blocks,threads>>>(recvcnt, recvmap.data(), recvmap.size(), recv_buffer, recv_buffer_size, dfield);

  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error launching GPU packing kernel")+ cudaGetErrorString(err);
    throw_Exception(msg);
  }

  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    std::string msg = std::string("Error synchronizing device")+ cudaGetErrorString(err);
    throw_Exception(msg);
  }
}


#define ATLAS_REPEAT_MACRO(n, m, p) ATLAS_REPEAT_MACRO_ ## n(m, p)
// expands to m(0,p) m(1,p) ... m(n-1,p)

#define ATLAS_REPEAT_MACRO_0(m, p)
#define ATLAS_REPEAT_MACRO_1(m, p) ATLAS_REPEAT_MACRO_0(m, p) m(0, p)
#define ATLAS_REPEAT_MACRO_2(m, p) ATLAS_REPEAT_MACRO_1(m, p) m(1, p)
#define ATLAS_REPEAT_MACRO_3(m, p) ATLAS_REPEAT_MACRO_2(m, p) m(2, p)
#define ATLAS_REPEAT_MACRO_4(m, p) ATLAS_REPEAT_MACRO_3(m, p) m(3, p)
#define ATLAS_REPEAT_MACRO_5(m, p) ATLAS_REPEAT_MACRO_4(m, p) m(4, p)
#define ATLAS_REPEAT_MACRO_6(m, p) ATLAS_REPEAT_MACRO_5(m, p) m(5, p)


#define EXPLICIT_TEMPLATE_INSTANTIATION( ParallelDim, RANK )         \
  template class halo_packer_cuda<ParallelDim, int,           RANK>; \
  template class halo_packer_cuda<ParallelDim, long,          RANK>; \
  template class halo_packer_cuda<ParallelDim, long unsigned, RANK>; \
  template class halo_packer_cuda<ParallelDim, float,         RANK>; \
  template class halo_packer_cuda<ParallelDim, double,        RANK>;

#define EXPLICIT_TEMPLATE_INSTANTIATION_RANK(RANK) \
  ATLAS_REPEAT_MACRO(RANK, EXPLICIT_TEMPLATE_INSTANTIATION, RANK)

EXPLICIT_TEMPLATE_INSTANTIATION_RANK(1)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(2)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(3)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(4)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(5)

} //namespace array
} //namespace atlas
