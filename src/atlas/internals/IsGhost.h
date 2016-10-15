/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_util_IsGhost_h
#define atlas_util_IsGhost_h

#include "atlas/mesh/Nodes.h"
#include "atlas/field/Field.h"
#include "atlas/internals/Bitflags.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace internals {

class IsGhost
{
public:
  IsGhost( const mesh::Nodes& nodes )
  {
//    part_   = array::ArrayView<int,1> (nodes.partition() );
//    ridx_   = array::IndexView<int,1> (nodes.remote_index() );
//    mypart_ = parallel::mpi::comm().rank();
  flags_ = array::ArrayView<int,1> (nodes.field("flags"));
  ghost_ = array::ArrayView<int,1> (nodes.ghost());
  }

  bool operator()(size_t idx) const
  {
//    if( part_[idx] != mypart_ ) return true;
//    if( ridx_[idx] != idx     ) return true;
//    return false;
    return Topology::check(flags_(idx),Topology::GHOST);
//    return ghost_(idx);
  }
private:
//  int mypart_;
//  array::ArrayView<int,1> part_;
//  array::IndexView<int,1> ridx_;
  array::ArrayView<int,1> flags_;
  array::ArrayView<int,1> ghost_;
};

} // namespace internals
} // namespace atlas

#endif
