/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_IsGhost_h
#define atlas_util_IsGhost_h

#include "atlas/mpi/mpi.h"
#include "atlas/FunctionSpace.h"
#include "atlas/Field.h"
#include "atlas/util/ArrayView.h"
#include "atlas/util/IndexView.h"

namespace atlas {
namespace util {

struct IsGhost
{
  IsGhost( FunctionSpace& nodes )
  {
    part_   = ArrayView<int,1> (nodes.field("partition") );
    ridx_   = IndexView<int,1> (nodes.field("remote_idx") );
    mypart_ = eckit::mpi::rank();
  }
  IsGhost( FunctionSpace& nodes, int mypart )
  {
    part_   = ArrayView<int,1> (nodes.field("partition") );
    ridx_   = IndexView<int,1> (nodes.field("remote_idx") );
    mypart_ = mypart;
  }

  bool operator()(int idx)
  {
    if( part_[idx] != mypart_ ) return true;
    if( ridx_[idx] != idx     ) return true;
    return false;
  }
  int mypart_;
  ArrayView<int,1> part_;
  IndexView<int,1> ridx_;
};

} // namespace util
} // namespace atlas

#endif
