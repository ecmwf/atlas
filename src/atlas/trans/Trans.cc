/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grids/ReducedGrid.h"
#include "atlas/ErrorHandling.h"
#include "atlas/trans/Trans.h"

namespace atlas {
namespace trans {

Trans::Trans(const grids::ReducedGrid& g)
{
  int nsmax = (2*g.nlat()-1)/2;
  ctor(g.nlat(),g.npts_per_lat().data(), nsmax);
}

Trans::Trans(const grids::ReducedGrid& g, const int nsmax )
{
  ctor(g.nlat(),g.npts_per_lat().data(), nsmax);
}

Trans::Trans( const std::vector<int>& npts_per_lat, const int nsmax )
{
  ctor(npts_per_lat.size(),npts_per_lat.data(), nsmax);
}

Trans::~Trans()
{
  ::trans_delete(&trans_);
}

void Trans::ctor(const int ndgl, const int nloen[], int nsmax)
{
  trans_ = ::new_trans();
  trans_.ndgl  = ndgl;
  trans_.nloen = new int[trans_.ndgl];
  std::copy(nloen,nloen+ndgl,trans_.nloen);
  trans_.nsmax = nsmax;
  if( nsmax == 0 )
    trans_.lgridonly = true;
  trans_setup(&trans_);
}


Trans* atlas__Trans__new (grids::ReducedGrid* grid)
{
  Trans* trans;
  ATLAS_ERROR_HANDLING(
    ASSERT( grid != NULL );
    trans = new Trans(*grid);
  );
  return trans;
}

void atlas__Trans__delete (Trans* trans)
{
  ATLAS_ERROR_HANDLING( delete trans );
}

int atlas__Trans__handle (Trans* trans)
{
  ATLAS_ERROR_HANDLING( return trans->handle() );
}

}
}
