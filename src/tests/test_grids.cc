/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <sstream>
#include <algorithm>
#include <iomanip>
#define BOOST_TEST_MODULE TestGrids
#include "ecbuild/boost_test_framework.h"

#include <eckit/memory/Factory.h>
#include <eckit/memory/Builder.h>
#include "atlas/atlas.h"
#include "atlas/ReducedGrid.h"

using namespace eckit;
using namespace atlas;

BOOST_AUTO_TEST_CASE( init ) { MPL::init(); }

BOOST_AUTO_TEST_CASE( test_factory )
{
  std::cout << Factory<ReducedGrid>::instance() << std::endl;
  SharedPtr<ReducedGrid> reducedgrid( Factory<ReducedGrid>::instance().get("reduced_gg.N80").create() );
  SharedPtr<Grid> grid( Factory<Grid>::instance().get("reduced_gg.N24").create(eckit::ValueParams()) );
  std::cout << "reducedgrid->nlat() = " << reducedgrid->nlat() << std::endl;
  std::cout << "grid->nPoints() = " << grid->nPoints() << std::endl;

  std::vector<std::string> keys = Factory<ReducedGrid>::instance().keys();
  for( int i=0; i< keys.size(); ++i )
  {
    std::cout << keys[i] << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( finalize ) { MPL::finalize(); }
