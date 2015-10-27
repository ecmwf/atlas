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
#include <iterator>     // std::iterator, std::input_iterator_tag

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
#include "atlas/util/IndexView.h"

// ------------------------------------------------------------------

using namespace atlas::mesh;

namespace atlas {
namespace test {

class Quadrilateral : public ElementType
{
public:
  virtual ~Quadrilateral() {}
  virtual size_t nb_nodes() const { return 4; }
  virtual size_t nb_edges() const { return 4; }
  virtual const std::string& name() const { static std::string s("Quadrilateral"); return s; }
};

class Triangle : public ElementType
{
public:
  virtual ~Triangle() {}
  virtual size_t nb_nodes() const { return 3; }
  virtual size_t nb_edges() const { return 3; }
  virtual const std::string& name() const { static std::string s("Triangle"); return s; }
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
  Nodes nodes(6);
  Elements elements(nodes);

  size_t triags[] = {
    1,5,3,
    1,5,2
  };
  elements.add( new Triangle(), 2, triags );

  size_t quads[] = {
    0,1,2,3
  };
  elements.add( new Quadrilateral(), 1, quads );

  {
    const Connectivity& connectivity = elements.node_connectivity();
    for( size_t e=0; e<elements.size(); ++e ) {
      eckit::Log::info() << e << std::endl;
      eckit::Log::info() << "  " << elements.name(e) << std::endl;
      eckit::Log::info() << "  nb_nodes = " << elements.nb_nodes(e) << std::endl;
      eckit::Log::info() << "  nb_edges = " << elements.nb_edges(e) << std::endl;
      eckit::Log::info() << "  nodes = [ ";
      for( size_t n=0; n<elements.nb_nodes(e); ++n ) {
        eckit::Log::info() << connectivity[e][n] << " ";
      }
      eckit::Log::info() << "]" << std::endl;
    }
  }

  eckit::Log::info() << std::endl;

  {
    for( size_t t=0; t<elements.nb_types(); ++t ) {
      const ElementTypeElements etype_elements(elements,t);
      eckit::Log::info() << "name = " << etype_elements.name() << std::endl;
      eckit::Log::info() << "nb_elements = " << etype_elements.size() << std::endl;
      const ElementTypeConnectivity& connectivity = etype_elements.node_connectivity();
      for( size_t e=0; e<elements.nb_elements(t); ++e ) {
        eckit::Log::info() << "  nodes = [ ";
        for( size_t n=0; n<etype_elements.nb_nodes(e); ++n ) {
          eckit::Log::info() << connectivity[e][n] << " ";
        }
        eckit::Log::info() << "]" << std::endl;
      }
    }
  }

}

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
