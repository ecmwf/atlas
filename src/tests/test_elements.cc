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
#include "atlas/Connectivity.h"

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

BOOST_AUTO_TEST_CASE( hybrid_elements )
{
  HybridElements hybrid_elements;

  idx_t triangle_nodes[] = {
    1,5,3,
    1,5,2
  };
  size_t triags = hybrid_elements.add( new Triangle(), 2, triangle_nodes );

  std::vector<idx_t> quad_nodes(4);
  quad_nodes[0] = 0;
  quad_nodes[1] = 1;
  quad_nodes[2] = 2;
  quad_nodes[3] = 3;
  size_t quads = hybrid_elements.add( new Quadrilateral(), 1, quad_nodes );


  {
    const HybridElements::Connectivity& connectivity = hybrid_elements.node_connectivity();
    for( size_t e=0; e<hybrid_elements.size(); ++e ) {
      eckit::Log::info() << e << std::endl;
      eckit::Log::info() << "  " << hybrid_elements.name(e) << std::endl;
      eckit::Log::info() << "  nb_nodes = " << hybrid_elements.nb_nodes(e) << std::endl;
      eckit::Log::info() << "  nb_edges = " << hybrid_elements.nb_edges(e) << std::endl;
      eckit::Log::info() << "  nodes = [ ";
      for( size_t n=0; n<hybrid_elements.nb_nodes(e); ++n ) {
        eckit::Log::info() << connectivity(e,n) << " ";
      }
      eckit::Log::info() << "]" << std::endl;
    }
    BOOST_CHECK_EQUAL( connectivity(0,0) , 1 );
    BOOST_CHECK_EQUAL( connectivity(0,1) , 5 );
    BOOST_CHECK_EQUAL( connectivity(0,2) , 3 );
    BOOST_CHECK_EQUAL( connectivity(1,0) , 1 );
    BOOST_CHECK_EQUAL( connectivity(1,1) , 5 );
    BOOST_CHECK_EQUAL( connectivity(1,2) , 2 );
    BOOST_CHECK_EQUAL( connectivity(2,0) , 0 );
    BOOST_CHECK_EQUAL( connectivity(2,1) , 1 );
    BOOST_CHECK_EQUAL( connectivity(2,2) , 2 );
    BOOST_CHECK_EQUAL( connectivity(2,3) , 3 );
  }

  eckit::Log::info() << std::endl;
  idx_t quad0[4] = {9,8,7,6};
  {
    for( size_t t=0; t<hybrid_elements.nb_types(); ++t ) {
      Elements& elements = hybrid_elements.elements(t);
      const Elements::Connectivity& block_connectivity = hybrid_elements.node_connectivity(t);
      if( t==0 ) {
        BOOST_CHECK_EQUAL( block_connectivity(0,0) , 1 );
        BOOST_CHECK_EQUAL( block_connectivity(0,1) , 5 );
        BOOST_CHECK_EQUAL( block_connectivity(0,2) , 3 );
        BOOST_CHECK_EQUAL( block_connectivity(1,0) , 1 );
        BOOST_CHECK_EQUAL( block_connectivity(1,1) , 5 );
        BOOST_CHECK_EQUAL( block_connectivity(1,2) , 2 );
      }
      if( t==1 ) {
        BOOST_CHECK_EQUAL( block_connectivity(0,0) , 0 );
        BOOST_CHECK_EQUAL( block_connectivity(0,1) , 1 );
        BOOST_CHECK_EQUAL( block_connectivity(0,2) , 2 );
        BOOST_CHECK_EQUAL( block_connectivity(0,3) , 3 );
      }
      const_cast<Elements::Connectivity&>(block_connectivity).set(0,quad0);
      elements.set_node_connectivity(0,quad0);

      if( t==0 ) {
        BOOST_CHECK_EQUAL( block_connectivity(0,0) , 9 );
        BOOST_CHECK_EQUAL( block_connectivity(0,1) , 8 );
        BOOST_CHECK_EQUAL( block_connectivity(0,2) , 7 );
        BOOST_CHECK_EQUAL( block_connectivity(1,0) , 1 );
        BOOST_CHECK_EQUAL( block_connectivity(1,1) , 5 );
        BOOST_CHECK_EQUAL( block_connectivity(1,2) , 2 );
      }
      if( t==1 ) {
        BOOST_CHECK_EQUAL( block_connectivity(0,0) , 9 );
        BOOST_CHECK_EQUAL( block_connectivity(0,1) , 8 );
        BOOST_CHECK_EQUAL( block_connectivity(0,2) , 7 );
        BOOST_CHECK_EQUAL( block_connectivity(0,3) , 6 );
      }

      eckit::Log::info() << "name = " << elements.name() << std::endl;
      eckit::Log::info() << "nb_elements = " << elements.size() << std::endl;
      const Elements::Connectivity& connectivity = elements.node_connectivity();
      for( size_t e=0; e<elements.size(); ++e ) {
        eckit::Log::info() << "  nodes = [ ";
        for( size_t n=0; n<elements.nb_nodes(); ++n ) {
          eckit::Log::info() << connectivity(e,n) << " ";
        }
        eckit::Log::info() << "]" << std::endl;
      }
    }
  }

}

BOOST_AUTO_TEST_CASE( elements )
{
  eckit::Log::info() << "\nelements" << std::endl;

  idx_t triangle_nodes[] = {
    1,5,3,
    1,5,2
  };
  idx_t triag1[3] = {9,8,7};

  Elements elements( new Triangle(), 2, triangle_nodes );
  elements.set_node_connectivity(0,triag1);
  eckit::Log::info() << "name = " << elements.name() << std::endl;
  eckit::Log::info() << "nb_elements = " << elements.size() << std::endl;
  const Elements::Connectivity& connectivity = elements.node_connectivity();
  for( size_t e=0; e<elements.size(); ++e ) {
    eckit::Log::info() << "  nodes = [ ";
    for( size_t n=0; n<elements.nb_nodes(); ++n ) {
       eckit::Log::info() << connectivity(e,n) << " ";
    }
    eckit::Log::info() << "]" << std::endl;
  }

  // This will create a new structure with a copy of elements
  HybridElements hybrid_elements;
  hybrid_elements.add(elements);
  {
    for( size_t t=0; t<hybrid_elements.nb_types(); ++t ) {
      Elements& elements = hybrid_elements.elements(t);
      elements.set_node_connectivity(0,triag1);
      eckit::Log::info() << "name = " << elements.name() << std::endl;
      eckit::Log::info() << "nb_elements = " << elements.size() << std::endl;
      const Elements::Connectivity& connectivity = elements.node_connectivity();
      for( size_t e=0; e<elements.size(); ++e ) {
        eckit::Log::info() << "  nodes = [ ";
        for( size_t n=0; n<elements.nb_nodes(); ++n ) {
          eckit::Log::info() << connectivity(e,n) << " ";
        }
        eckit::Log::info() << "]" << std::endl;
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( hybrid_connectivity )
{
  eckit::Log::info() << "\nhybrid_connectivity" << std::endl;

  idx_t triangle_nodes[] = {
    1,5,3,
    1,5,2
  };
  MultiBlockConnectivity hybrid_connectivity;
  hybrid_connectivity.add(2,3,triangle_nodes);
  for( size_t e=0; e<hybrid_connectivity.rows(); ++e )
  {
    eckit::Log::info() << "  cols = [ ";
    for( size_t n=0; n<hybrid_connectivity.cols(e); ++n ) {
      eckit::Log::info() << hybrid_connectivity(e,n) << " ";
    }
    eckit::Log::info() << "]" << std::endl;
  }


  idx_t quad_nodes[] = { 0,1,2,3,
                         4,5,6,7,
                         8,9,10,11 };
  hybrid_connectivity.add(3,4, quad_nodes);


  for( size_t e=0; e<hybrid_connectivity.rows(); ++e )
  {
    eckit::Log::info() << "  cols = [ ";
    for( size_t n=0; n<hybrid_connectivity.cols(e); ++n ) {
      eckit::Log::info() << hybrid_connectivity(e,n) << " ";
    }
    eckit::Log::info() << "]" << std::endl;
  }

  for( size_t b=0; b<hybrid_connectivity.blocks(); ++b )
  {
    const BlockConnectivity& block = hybrid_connectivity.block_connectivity(b);
    for( size_t r=0; r<block.rows(); ++r )
    {
      eckit::Log::info() << "  cols = [ ";
      for( size_t c=0; c<block.cols(); ++c ) {
        eckit::Log::info() << block(r,c) << " ";
      }
      eckit::Log::info() << "]" << std::endl;
    }
  }


}


BOOST_AUTO_TEST_CASE( block_connectivity )
{
  eckit::Log::info() << "\nblock_connectivity" << std::endl;

  idx_t triangle_nodes[] = {
    1,5,3,
    1,5,2
  };
  BlockConnectivity block;
  block.add(2,3,triangle_nodes);
  block.add(2,3,triangle_nodes);
  for( size_t r=0; r<block.rows(); ++r )
  {
    eckit::Log::info() << "  cols = [ ";
    for( size_t c=0; c<block.cols(); ++c ) {
      eckit::Log::info() << block(r,c) << " ";
    }
    eckit::Log::info() << "]" << std::endl;
  }

}



BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
