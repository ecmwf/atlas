/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <eckit/log/Log.h>
#include "atlas/grids/grids.h"


namespace atlas {
namespace grids {

template<typename CONCRETE>
void load_grid()
{
  eckit::ConcreteBuilderT1<Grid,CONCRETE> builder("tmp");
}

void load()
{
  eckit::Log::debug() << "Loading library [atlas::grids]" << std::endl;

  // We have to touch all classes we want to register for static linking.

  load_grid<ReducedGrid>();
  load_grid<GaussianGrid>();
  load_grid<ReducedGaussianGrid>();
  load_grid<LonLatGrid>();
  load_grid<ReducedLonLatGrid>();
  load_grid<RotatedLatLon>();
  load_grid<Unstructured>();
  load_grid<PolarStereoGraphic>();

  load_grid<rgg::N16>();
  load_grid<rgg::N24>();
  load_grid<rgg::N32>();
  load_grid<rgg::N48>();
  load_grid<rgg::N64>();
  load_grid<rgg::N80>();
  load_grid<rgg::N96>();
  load_grid<rgg::N128>();
  load_grid<rgg::N160>();
  load_grid<rgg::N200>();
  load_grid<rgg::N256>();
  load_grid<rgg::N320>();
  load_grid<rgg::N400>();
  load_grid<rgg::N512>();
  load_grid<rgg::N576>();
  load_grid<rgg::N640>();
  load_grid<rgg::N800>();
  load_grid<rgg::N1024>();
  load_grid<rgg::N1280>();
  load_grid<rgg::N1600>();
  load_grid<rgg::N2000>();
  load_grid<rgg::N4000>();
  load_grid<rgg::N8000>();
}

} // namespace grids
} // namespace atlas

extern "C"
{
  void atlas__grids__load()
  {
    atlas::grids::load();
  }
}
