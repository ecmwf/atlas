/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/mpl/MPL.hpp"
#include "atlas/atlas_config.h"
#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/FunctionSpace.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"
#include "atlas/mesh/Parameters.hpp"

using namespace atlas;
using namespace atlas::meshgen;

namespace atlas {

struct IsGhost
{
  IsGhost( FunctionSpace& nodes )
  {
    part    = ArrayView<int,1> (nodes.field("partition") );
    loc_idx = IndexView<int,1> (nodes.field("remote_idx") );
    mypart  = MPL::rank();
  }
  bool operator()(int idx)
  {
    if( part   [idx] != mypart ) return true;
    if( loc_idx[idx] != idx    ) return true;
    return false;
  }
  int mypart;
  ArrayView<int,1> part;
  IndexView<int,1> loc_idx;
};

} // namespace atlas
