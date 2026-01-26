/*
 * (C) Copyright 2025- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "Packer.h"

#include "pack_index.h"
#include "adjoint_unpack_index.h"
#include "zero_index.h"

namespace atlas::parallel::detail {

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HostPacker<ParallelDim,DATA_TYPE,RANK>::pack(const int sendcnt, array::SVector<int> const& sendmap,
                      const array::ArrayView<DATA_TYPE, RANK>& array, DATA_TYPE* send_buffer, int /*send_buffer_size*/) {
    idx_t ibuf = 0;
    for (int node_cnt = 0; node_cnt < sendcnt; ++node_cnt) {
        const idx_t node_idx = sendmap[node_cnt];
        pack_index<ParallelDim, RANK, 0>::apply(ibuf, node_idx, array, send_buffer);
    }
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HostPacker<ParallelDim,DATA_TYPE,RANK>::unpack(const int recvcnt, array::SVector<int> const& recvmap, const DATA_TYPE* recv_buffer,
                        int /*recv_buffer_size*/, array::ArrayView<DATA_TYPE, RANK>& array) {
    idx_t ibuf = 0;
    for (int node_cnt = 0; node_cnt < recvcnt; ++node_cnt) {
        const idx_t node_idx = recvmap[node_cnt];
        unpack_index<ParallelDim, RANK, 0>::apply(ibuf, node_idx, recv_buffer, array);
    }
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HostAdjointPacker<ParallelDim,DATA_TYPE,RANK>::unpack(const int recvcnt, array::SVector<int> const& recvmap, const DATA_TYPE* recv_buffer,
                        int /*recv_buffer_size*/, array::ArrayView<DATA_TYPE, RANK>& array) {
    idx_t ibuf = 0;
    for (int node_cnt = 0; node_cnt < recvcnt; ++node_cnt) {
        const idx_t node_idx = recvmap[node_cnt];
        adjoint_unpack_index<ParallelDim, RANK, 0>::apply(ibuf, node_idx, recv_buffer, array);
    }
}

template <int ParallelDim, typename DATA_TYPE, int RANK>
void HostZeroer<ParallelDim,DATA_TYPE,RANK>::zero(const int sendcnt, array::SVector<int> const& sendmap, array::ArrayView<DATA_TYPE, RANK>& array,
                       DATA_TYPE* recv_buffer, int recv_buffer_size) {
    idx_t ibuf = 0;
    for (int node_cnt = 0; node_cnt < sendcnt; ++node_cnt) {
        const idx_t node_idx = sendmap[node_cnt];
        zero_index<ParallelDim, RANK, 0>::apply(ibuf, node_idx, array, recv_buffer);
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


#define EXPLICIT_TEMPLATE_INSTANTIATION( ParallelDim, RANK )           \
  template struct HostPacker<ParallelDim, int,           RANK>;        \
  template struct HostPacker<ParallelDim, long,          RANK>;        \
  template struct HostPacker<ParallelDim, long unsigned, RANK>;        \
  template struct HostPacker<ParallelDim, float,         RANK>;        \
  template struct HostPacker<ParallelDim, double,        RANK>;        \
  template struct HostAdjointPacker<ParallelDim, int,           RANK>; \
  template struct HostAdjointPacker<ParallelDim, long,          RANK>; \
  template struct HostAdjointPacker<ParallelDim, long unsigned, RANK>; \
  template struct HostAdjointPacker<ParallelDim, float,         RANK>; \
  template struct HostAdjointPacker<ParallelDim, double,        RANK>; \
  template struct HostZeroer<ParallelDim, int,           RANK>;        \
  template struct HostZeroer<ParallelDim, long,          RANK>;        \
  template struct HostZeroer<ParallelDim, long unsigned, RANK>;        \
  template struct HostZeroer<ParallelDim, float,         RANK>;        \
  template struct HostZeroer<ParallelDim, double,        RANK>;

#define EXPLICIT_TEMPLATE_INSTANTIATION_RANK(RANK) \
  ATLAS_REPEAT_MACRO(RANK, EXPLICIT_TEMPLATE_INSTANTIATION, RANK)

EXPLICIT_TEMPLATE_INSTANTIATION_RANK(1)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(2)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(3)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(4)
EXPLICIT_TEMPLATE_INSTANTIATION_RANK(5)

} //namespace atlas::array::detail
