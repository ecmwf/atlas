/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/array/ArrayView.h"
#include "atlas/array/ArrayViewDefs.h"
#include "atlas/array/SVector.h"

namespace atlas {
namespace parallel {

template <int ParallelDim, typename DATA_TYPE, int RANK>
struct halo_packer_cuda {
    static void pack( const int sendcnt, array::SVector<int> const& sendmap,
                      const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield,
                      const array::ArrayView<DATA_TYPE, RANK>& dfield, array::SVector<DATA_TYPE>& send_buffer );

    static void unpack( const int sendcnt, array::SVector<int> const& recvmap,
                        const array::SVector<DATA_TYPE>& recv_buffer,
                        const array::ArrayView<DATA_TYPE, RANK, array::Intent::ReadOnly>& hfield,
                        array::ArrayView<DATA_TYPE, RANK>& dhfield );
};

}  // namespace parallel
}  // namespace atlas
