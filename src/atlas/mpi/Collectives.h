/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef ATLAS_MPI_COLLECTIVES_h
#define ATLAS_MPI_COLLECTIVES_h

#include "atlas/mpi/mpi.h"
#include "atlas/util/ArrayView.h"

namespace atlas {
namespace mpi {

/// @brief Buffer<DATA_TYPE,SHAPE>
///
/// Class that keeps allocation of a MPI buffer including
/// counts and displacements, but with added index operator[]
/// that returns an ArrayView<DATA_TYPE,SHAPE> of the part
/// of the buffer for a processor index.
template <typename DATA_TYPE,int SHAPE> struct Buffer;

// ----------------------------------------------------------------------------------

template <typename DATA_TYPE,int SHAPE>
struct Buffer : eckit::mpi::Buffer<DATA_TYPE>
{
};

template <typename DATA_TYPE>
struct Buffer<DATA_TYPE,1> : public eckit::mpi::Buffer<DATA_TYPE>
{
  ArrayView<DATA_TYPE,1> operator[](int p)
  {
    return ArrayView<DATA_TYPE,1> ( eckit::mpi::Buffer<DATA_TYPE>::buf.data()+eckit::mpi::Buffer<DATA_TYPE>::displs[p],
                                    make_shape( eckit::mpi::Buffer<DATA_TYPE>::counts[p] ).data() );
  }
};

// ----------------------------------------------------------------------------------

} // namespace mpi
} // namepsace atlas

#endif // ATLAS_MPI_COLLECTIVES_h
