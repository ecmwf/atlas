/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <cmath>
#include <string>
#include <iterator>

#include "atlas/internals/atlas_config.h"

#define BOOST_TEST_MODULE atlas_test_elements
#include "ecbuild/boost_test_framework.h"

#include "eckit/memory/ScopedPtr.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/atlas.h"
#include "atlas/runtime/Log.h"
#include "atlas/mesh/ElementType.h"
#include "atlas/mesh/Elements.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/field/Field.h"
#include "atlas/mesh/Connectivity.h"

#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/generators/Structured.h"
#include "atlas/grid/Grid.h"

#include "tests/AtlasFixture.h"


// ------------------------------------------------------------------

using namespace atlas::mesh;
using namespace atlas::mesh::temporary;

namespace atlas {
namespace test {


// ===================================================================
//                               BEGIN TESTS
// ===================================================================


BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( hybrid_elements )
{
  HybridElements hybrid_elements;

  idx_t triangle_nodes[] = {
    1,5,3,
    1,5,2
  };

  size_t triags_type_idx = hybrid_elements.add( new Triangle(), 2, triangle_nodes );

  BOOST_CHECK_EQUAL( triags_type_idx , 0 );

  hybrid_elements.add(field::Field::create<double>("surface",array::make_shape(hybrid_elements.size())));

  std::vector<idx_t> quad_nodes(4);
  quad_nodes[0] = 0;
  quad_nodes[1] = 1;
  quad_nodes[2] = 2;
  quad_nodes[3] = 3;

  size_t quads_type_idx = hybrid_elements.add( new Quadrilateral(), 1, quad_nodes );

  BOOST_CHECK_EQUAL( quads_type_idx , 1 );

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
      const Elements::Connectivity& block_connectivity = elements.node_connectivity();
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
      elements.node_connectivity().set(0,quad0);
      //const_cast<Elements::Connectivity&>(block_connectivity).set(0,quad0);
      //elements.set_node_connectivity(0,quad0);

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

  size_t nb_elements = 3;
  BOOST_CHECK_EQUAL( hybrid_elements.size(),                  nb_elements);
  BOOST_CHECK_EQUAL( hybrid_elements.global_index().size(),   nb_elements);
  BOOST_CHECK_EQUAL( hybrid_elements.partition().size(),      nb_elements);
  BOOST_CHECK_EQUAL( hybrid_elements.remote_index().size(),   nb_elements);
  BOOST_CHECK_EQUAL( hybrid_elements.field("surface").size(), nb_elements);
  BOOST_CHECK_EQUAL( hybrid_elements.elements(0).begin(), 0);
  BOOST_CHECK_EQUAL( hybrid_elements.elements(0).end(),   2);
  BOOST_CHECK_EQUAL( hybrid_elements.elements(1).begin(), 2);
  BOOST_CHECK_EQUAL( hybrid_elements.elements(1).end(),   3);
//  BOOST_CHECK_EQUAL( hybrid_elements.elements(0).type_idx(), 0 );
//  BOOST_CHECK_EQUAL( hybrid_elements.elements(1).type_idx(), 1 );
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

  BOOST_CHECK_EQUAL( elements.begin(), 0);
  BOOST_CHECK_EQUAL( elements.end(),   2);
//  BOOST_CHECK_EQUAL( elements.type_idx(), 0);

  elements.node_connectivity().set(0,triag1);
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
      elements.node_connectivity().set(0,triag1);
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
    const BlockConnectivity& block = hybrid_connectivity.block(b);
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

BOOST_AUTO_TEST_CASE( zero_elements )
{
  HybridElements hybrid_elements;
  idx_t *nodes = 0;

  hybrid_elements.add(new Triangle(), 0, nodes );
  hybrid_elements.add(new Quadrilateral(), 0, nodes );

  BOOST_CHECK_EQUAL( hybrid_elements.size(), 0 );
  BOOST_CHECK_EQUAL( hybrid_elements.nb_types(), 2 );
  BOOST_CHECK_EQUAL( hybrid_elements.elements(0).size(), 0 );
  BOOST_CHECK_EQUAL( hybrid_elements.elements(1).size(), 0 );
}

BOOST_AUTO_TEST_CASE( irregularconnectivity_insert )
{
  IrregularConnectivity connectivity;
  idx_t c[] = {1, 2, 3, 4,
               5, 6, 7, 8,
               9,10,11,12};
  connectivity.add(3,4,c);
  idx_t c2[] = {13,14,15,
                16,17,18};
  connectivity.insert(1,2,3,c2);
  connectivity.insert(2,1,5);

  size_t iregular_c[] = {2, 3, 4, 1};
  connectivity.insert(5,4,iregular_c);

  idx_t values[]={ 1, 2, 3, 4, 13, 14, 15, -1, -1, -1, -1, -1, 16, 17, 18, 5, 6, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 9, 10, 11, 12 };
  idx_t counts[]={ 4, 3, 5, 3, 4, 2, 3, 4, 1, 4 };

  size_t n(0);
  size_t r(0);
  for( size_t jrow=0; jrow<connectivity.rows(); ++jrow )
  {
    BOOST_CHECK_EQUAL(connectivity.cols(jrow), counts[r++]);
    for( size_t jcol=0; jcol<connectivity.cols(jrow); ++jcol )
    {
      BOOST_CHECK_EQUAL(connectivity(jrow,jcol), values[n++]);
    }
  }
}

BOOST_AUTO_TEST_CASE( multiblockconnectivity_insert )
{
  MultiBlockConnectivity connectivity;
  idx_t c1[] = {1, 2, 3, 4,
                5, 6, 7, 8,
                9,10,11,12};
  connectivity.add(3,4,c1);
  idx_t c2[] = {13, 14, 15,
                16, 17, 18};
  connectivity.add(2,3,c2);

  idx_t c1i[] = {19,20,21,22};
  idx_t c2i[] = {23,24,25,
                 26,27,28};

  BOOST_CHECK_EQUAL( connectivity.block(0).rows() , 3 );
  BOOST_CHECK_EQUAL( connectivity.block(1).rows() , 2 );
  BOOST_CHECK_EQUAL( connectivity.block(0).cols() , 4 );
  BOOST_CHECK_EQUAL( connectivity.block(1).cols() , 3 );

  std::cout << "block 0" << std::endl;
  for( size_t jrow=0; jrow<connectivity.block(0).rows(); ++jrow )
  {
    for( size_t jcol=0; jcol<connectivity.block(0).cols(); ++jcol )
    {
      std::cout << connectivity.block(0)(jrow,jcol) << " ";
    }
    std::cout << std::endl;
  }

  connectivity.insert(1,1,4,c1i);
  connectivity.insert(5,2,3,c2i);

  BOOST_CHECK_EQUAL( connectivity.block(0).rows() , 4 );
  BOOST_CHECK_EQUAL( connectivity.block(1).rows() , 4 );

  std::cout << "\nfull\n";
  for( size_t jrow=0; jrow<connectivity.rows(); ++jrow )
  {
    for( size_t jcol=0; jcol<connectivity.cols(jrow); ++jcol )
    {
      std::cout << connectivity(jrow,jcol) << " ";
    }
    std::cout << std::endl;
  }
  idx_t values[]={ 1, 2, 3, 4, 19, 20, 21, 22, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 23, 24, 25, 26, 27, 28, 16, 17, 18 };
  idx_t counts[]={ 4, 4, 4, 4, 3, 3, 3, 3 };

  size_t n(0);
  size_t r(0);
  for( size_t jrow=0; jrow<connectivity.rows(); ++jrow )
  {
    BOOST_CHECK_EQUAL(connectivity.cols(jrow), counts[r++]);
    for( size_t jcol=0; jcol<connectivity.cols(jrow); ++jcol )
    {
      BOOST_CHECK_EQUAL(connectivity(jrow,jcol), values[n++]);
    }
  }
}

BOOST_AUTO_TEST_CASE( cells_insert )
{
  Log::info() << "\n\n\ncells_insert \n============ \n\n" << std::endl;
  HybridElements cells;
  idx_t c1[] = {1, 2, 3, 4,
                5, 6, 7, 8,
                9,10,11,12};
  cells.add(new Quadrilateral(), 3, c1);
  idx_t c2[] = {13, 14, 15,
                16, 17, 18};
  cells.add(new Triangle(), 2, c2);

  BOOST_CHECK_EQUAL( cells.elements(0).size() , 3 );
  BOOST_CHECK_EQUAL( cells.elements(1).size() , 2 );
  BOOST_CHECK_EQUAL( cells.size(), 5 );
  BOOST_CHECK_EQUAL( cells.elements(0).node_connectivity().rows(), 3 );
  BOOST_CHECK_EQUAL( cells.elements(1).node_connectivity().rows(), 2 );

  Log::info() << "Update elements(0)" << std::endl;
  size_t pos0 = cells.elements(0).add(3);


  Log::info() << "Update elements(1)" << std::endl;
  size_t pos1 = cells.elements(1).add(2);

  BOOST_CHECK_EQUAL( pos0, 3 );
  BOOST_CHECK_EQUAL( pos1, 2 );
  BOOST_CHECK_EQUAL( cells.elements(0).size() , 6 );
  BOOST_CHECK_EQUAL( cells.elements(1).size() , 4 );
  BOOST_CHECK_EQUAL( cells.size(), 10 );
  BOOST_CHECK_EQUAL( cells.node_connectivity().block(0).rows(), 6 );
  BOOST_CHECK_EQUAL( cells.elements(0).node_connectivity().rows(), 6 );
  BOOST_CHECK_EQUAL( cells.elements(1).node_connectivity().rows(), 4 );

  const BlockConnectivity& conn1 = cells.elements(0).node_connectivity();
  const BlockConnectivity& conn2 = cells.elements(1).node_connectivity();

  std::cout << "\nconn1\n";
  for( size_t jrow=0; jrow<conn1.rows(); ++jrow )
  {
    for( size_t jcol=0; jcol<conn1.cols(); ++jcol )
    {
      std::cout << conn1(jrow,jcol) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "\nconn2\n";
  for( size_t jrow=0; jrow<conn2.rows(); ++jrow )
  {
    for( size_t jcol=0; jcol<conn2.cols(); ++jcol )
    {
      std::cout << conn2(jrow,jcol) << " ";
    }
    std::cout << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( cells_add_add )
{
  HybridElements cells;

  cells.add(new Quadrilateral(), 3);
  cells.add(new Triangle(),      2);

  HybridElements::Connectivity& conn = cells.node_connectivity();

  int nodes[] = {0,1,2,3};

  conn.set(0,nodes);

  BOOST_CHECK_EQUAL( conn(0,0), 0 );
  BOOST_CHECK_EQUAL( conn(0,1), 1 );
  BOOST_CHECK_EQUAL( conn(0,2), 2 );
  BOOST_CHECK_EQUAL( conn(0,3), 3 );

  BOOST_CHECK_EQUAL( cells.elements(0).node_connectivity()(0,0), 0 );
  BOOST_CHECK_EQUAL( cells.elements(0).node_connectivity()(0,1), 1 );
  BOOST_CHECK_EQUAL( cells.elements(0).node_connectivity()(0,2), 2 );
  BOOST_CHECK_EQUAL( cells.elements(0).node_connectivity()(0,3), 3 );

}


} // namespace test
} // namespace atlas
