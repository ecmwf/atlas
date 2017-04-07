/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestRotation
#include "ecbuild/boost_test_framework.h"


#include "atlas/library/Library.h"
#include "atlas/grid/detail/projection/Rotation.h"
#include "atlas/util/Config.h"
#include "atlas/runtime/Log.h"
#include "tests/AtlasFixture.h"

using atlas::grid::projection::Rotated;
using atlas::util::Config;

namespace atlas {
namespace test {

constexpr double eps() { return 1.e-5; }

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_rotation )
{
  Config config;
  config.set("north_pole", std::vector<double>{-176,40} );
  Rotated rotation(config);
  
  Log::info() << rotation << std::endl;
  
  PointLonLat p;
  
  p = {0.,90};
  rotation.rotate(p.data());
  BOOST_CHECK_CLOSE( p.lon(), -176. , eps() );
  BOOST_CHECK_CLOSE( p.lat(),   40. , eps() );
  rotation.unrotate(p.data());
  // p.lon() could be any value... singularity
  BOOST_CHECK_CLOSE( p.lat(),   90. , eps() );
 
  p = {0,0};
  rotation.rotate(p.data());
  BOOST_CHECK_CLOSE( p.lon(), -176. , eps() );
  BOOST_CHECK_CLOSE( p.lat(),  -50. , eps() );
  rotation.unrotate(p.data());
  BOOST_CHECK_SMALL( p.lon() , eps() );
  BOOST_CHECK_SMALL( p.lat() , eps() );

  p = {-180,45};
  rotation.rotate(p.data());
  BOOST_CHECK_CLOSE( p.lon(),  -176. , eps() );
  BOOST_CHECK_CLOSE( p.lat(),    85. , eps() );
  rotation.unrotate(p.data());
  BOOST_CHECK_CLOSE( p.lon(), -180. , eps() );
  BOOST_CHECK_CLOSE( p.lat(),   45. , eps() );
}


} // namespace test
} // namespace atlas
