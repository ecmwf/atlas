/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestMetadata
#include "ecbuild/boost_test_framework.h"

#include "atlas/atlas.h"
#include "atlas/internals/Debug.h"
#include "atlas/array/ArrayView.h"
#include "atlas/util/Metadata.h"

using namespace eckit;
using namespace atlas::util;

namespace atlas {
namespace test {

struct AtlasFixture {
    AtlasFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                 boost::unit_test::framework::master_test_suite().argv); }
    ~AtlasFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_broadcast_to_self )
{
  Metadata metadata;
  if( eckit::mpi::rank() == 0 )
  {
    metadata.set("paramID",128);
  }
  
  // broadcast
  metadata.broadcast();
  
  BOOST_CHECK( metadata.has("paramID") );
  if( metadata.has("paramID") )
    BOOST_CHECK_EQUAL( metadata.get<int>("paramID"), 128 );
  
}

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( test_broadcast_to_other )
{
  size_t root = 0;
  Metadata global;
  if( eckit::mpi::rank() == root )
  {
    global.set("paramID",128);
  }
  
  Metadata local;
  
  // broadcast
  global.broadcast(local);
  
  BOOST_CHECK( local.has("paramID") );
  if( local.has("paramID") )
    BOOST_CHECK_EQUAL( local.get<int>("paramID"), 128 );
  
  if( eckit::mpi::rank() != root )
    BOOST_CHECK( ! global.has("paramID") );
}

// -----------------------------------------------------------------------------

} // namespace test
} // namespace atlas
