/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "atlas/array/MakeView.h"
#include "atlas/parallel/mpi/mpi.h"

namespace atlas {
namespace internals {

class IsGhost
{
public:
  IsGhost( const mesh::Nodes& nodes ) :
    flags_( array::make_view<int,1> (nodes.field("flags")) ),
    ghost_( array::make_view<int,1> (nodes.ghost()) )
  {
  }

  bool operator()(size_t idx) const
  {
    return Topology::check(flags_(idx),Topology::GHOST);
  }
private:
  array::ArrayView<int,1> flags_;
  array::ArrayView<int,1> ghost_;
};

} // namespace internals
} // namespace atlas

#endif
