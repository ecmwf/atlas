/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>

#include "atlas/atlas_config.h"

#define BOOST_TEST_MODULE test_elements
#include "ecbuild/boost_test_framework.h"

#include "eckit/memory/ScopedPtr.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/atlas.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/Field.h"

// ------------------------------------------------------------------

using namespace atlas::mesh;

namespace atlas {
namespace test {

class MyElementType : public ElementType
{
  virtual ~MyElementType() {}
  virtual size_t nb_nodes() const { return 2; }
};

// ===================================================================
//                               BEGIN TESTS
// ===================================================================

struct GlobalFixture {
    GlobalFixture()  { atlas_init(boost::unit_test::framework::master_test_suite().argc,
                                  boost::unit_test::framework::master_test_suite().argv); }
    ~GlobalFixture() { atlas_finalize(); }
};

BOOST_GLOBAL_FIXTURE( GlobalFixture )

BOOST_AUTO_TEST_SUITE( test_elements )

BOOST_AUTO_TEST_CASE( simple )
{
  Nodes nodes(4);
  Elements elements(nodes);
  elements.add( new MyElementType(), 4 );

  size_t connectivity[] = {0,1,1,2};
  elements.add( new MyElementType(), 2, connectivity );
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
