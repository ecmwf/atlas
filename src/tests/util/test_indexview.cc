/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestIndexView
#include "ecbuild/boost_test_framework.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/IndexView.h"
#include "atlas/array/MakeView.h"
#include "tests/AtlasFixture.h"

#ifdef ATLAS_HAVE_FORTRAN
#define IN_FORTRAN
#else
#define IN_FORTRAN -1
#endif

using namespace atlas::array;

namespace atlas {
namespace test {

template< typename Iterator >
std::string pos(Iterator& it)
{
  std::stringstream ss;
  for(size_t i=0; i < it.pos().size(); ++i)  ss << it.pos()[i] << " ";
  return ss.str();
}

BOOST_GLOBAL_FIXTURE( AtlasFixture );
/*
BOOST_AUTO_TEST_CASE( test_array )
{
  array::ArrayT<int> _array (3,1,4);
  array::ArrayView<int,3> array(_array);
  BOOST_CHECK_EQUAL( array.shape(0) , 3 );
  BOOST_CHECK_EQUAL( array.shape(1) , 1 );
  BOOST_CHECK_EQUAL( array.shape(2) , 4 );
  BOOST_CHECK_EQUAL( array.size(), 12 );
  BOOST_CHECK_EQUAL( array.stride(0), 4 );
  BOOST_CHECK_EQUAL( array.stride(1), 4 );
  BOOST_CHECK_EQUAL( array.stride(2), 1 );
  for( size_t j=0; j<array.size(); ++j ) {
    array.data()[j] = j;
  }
  BOOST_CHECK_EQUAL( array(0,0,0) , 0 );
  BOOST_CHECK_EQUAL( array(0,0,1) , 1 );
  BOOST_CHECK_EQUAL( array(0,0,2) , 2 );
  BOOST_CHECK_EQUAL( array(0,0,3) , 3 );
  BOOST_CHECK_EQUAL( array(1,0,0) , 4 );
  BOOST_CHECK_EQUAL( array(1,0,1) , 5 );
  BOOST_CHECK_EQUAL( array(1,0,2) , 6 );
  BOOST_CHECK_EQUAL( array(1,0,3) , 7 );

#ifdef ATLAS_INDEXVIEW_BOUNDS_CHECKING
  BOOST_CHECK_THROW( array(0,1,0) , eckit::OutOfRange );  // j index out of range
  BOOST_CHECK_THROW( array(1,2,0,3), eckit::OutOfRange ); // rank out of range
#endif
}


BOOST_AUTO_TEST_CASE( test_arrayview_iterator )
{
  array::ArrayT<int> array(5,4,2);
  size_t strides[2] = {8,1};
  size_t extents[2] = {5,2};
  array::ArrayView<int>       aview(array.data(),2,extents,strides);
  array::ArrayView<int> const const_aview(array);

  std::cout << "aview.size() = " << aview.size() << std::endl;
  std::cout << "const_.size() = " << const_aview.size() << std::endl;

  array::ArrayView<int>::iterator it;
  array::ArrayView<int>::const_iterator const_it;

  int i(0);
  for(it = aview.begin(); it!=aview.end(); ++it, ++i)
  {
    std::cout << "set at pos " << test::pos(it) << " : " << i << std::endl;
    *it = i;
  }

  for( const_it = const_aview.begin(); const_it!=const_aview.end(); ++const_it)
  {
    std::cout << "read at pos " << test::pos(const_it) << " : " << *const_it << std::endl;
  }
}
*/

#ifdef ATLAS_HAVE_GRIDTOOLS_STORAGE
BOOST_AUTO_TEST_CASE( test_indexview_1d )
{
  //array::ArrayT<int> array( 10 );

  Array *array = Array::create<double>(10);
  ArrayView<double,1> aview = make_host_view<double, 1>(*array);

  aview(0) = 1 IN_FORTRAN;
  BOOST_CHECK_EQUAL( aview(0),       1 IN_FORTRAN );


 /*
  array::ArrayView<int,1>       aview(array);
  IndexView<int,1>       iview(array);
  const IndexView<int,1> const_iview(array);

  aview(0) = 1 IN_FORTRAN;
  BOOST_CHECK_EQUAL( aview(0),       1 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(0), 0 );
  BOOST_CHECK_EQUAL( int(iview(0)),  0 );

  iview(2) = 2;
  BOOST_CHECK_EQUAL( aview(2),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(2), 2 );
  BOOST_CHECK_EQUAL( int(iview(2)),  2 );

  iview(3) = iview(2);
  BOOST_CHECK_EQUAL( aview(3),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(3), 2 );
  BOOST_CHECK_EQUAL( int(iview(3)),  2 );

  iview(3) = iview(2) + 1;
  BOOST_CHECK_EQUAL( aview(3),       4 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(3), 3 );
  BOOST_CHECK_EQUAL( int(iview(3)) , 3 );

  iview(4) = iview(3);
  ++iview(4);
  BOOST_CHECK_EQUAL( aview(4),       5 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(4), 4 );
  BOOST_CHECK_EQUAL( int(iview(4)) , 4 );

  int val = iview(4);
  BOOST_CHECK_EQUAL( val, 4 );

  val = iview(4)+1;
  BOOST_CHECK_EQUAL( val, 5 );
*/
}
#else
BOOST_AUTO_TEST_CASE( test_indexview_1d )
{
  Array *array = Array::create<int>(10);

  ArrayView<int,1>       aview = make_view<int,1>(*array);
  IndexView<int,1>       iview = make_indexview<int,1>(*array);
  const IndexView<int,1> const_iview = make_indexview<int,1>(*array);

  aview(0) = 1 IN_FORTRAN;
  BOOST_CHECK_EQUAL( aview(0),       1 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(0), 0 );
  BOOST_CHECK_EQUAL( int(iview(0)),  0 );

  iview(2) = 2;
  BOOST_CHECK_EQUAL( aview(2),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(2), 2 );
  BOOST_CHECK_EQUAL( int(iview(2)),  2 );

  iview(3) = iview(2);
  BOOST_CHECK_EQUAL( aview(3),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(3), 2 );
  BOOST_CHECK_EQUAL( int(iview(3)),  2 );

  iview(3) = iview(2) + 1;
  BOOST_CHECK_EQUAL( aview(3),       4 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(3), 3 );
  BOOST_CHECK_EQUAL( int(iview(3)) , 3 );

  iview(4) = iview(3);
  ++iview(4);
  BOOST_CHECK_EQUAL( aview(4),       5 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(4), 4 );
  BOOST_CHECK_EQUAL( int(iview(4)) , 4 );

  int val = iview(4);
  BOOST_CHECK_EQUAL( val, 4 );

  val = iview(4)+1;
  BOOST_CHECK_EQUAL( val, 5 );

}
#endif

/*
BOOST_AUTO_TEST_CASE( test_indexview_2d )
{
  array::ArrayT<int> array( 5, 10 );

  array::ArrayView<int,2>       aview(array);
  IndexView<int,2>       iview(array);
  const IndexView<int,2> const_iview(array);

  int i=2;

  aview(i,0) = 1 IN_FORTRAN;
  BOOST_CHECK_EQUAL( aview(i,0),       1 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,0), 0 );
  BOOST_CHECK_EQUAL( int(iview(i,0)),  0 );

  iview(i,2) = 2;
  BOOST_CHECK_EQUAL( aview(i,2),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,2), 2 );
  BOOST_CHECK_EQUAL( int(iview(i,2)),  2 );

  iview(i,3) = iview(i,2);
  BOOST_CHECK_EQUAL( aview(i,3),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,3), 2 );
  BOOST_CHECK_EQUAL( int(iview(i,3)),  2 );

  iview(i,3) = iview(i,2) + 1;
  BOOST_CHECK_EQUAL( aview(i,3),       4 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,3), 3 );
  BOOST_CHECK_EQUAL( int(iview(i,3)) , 3 );

  iview(i,4) = iview(i,3);
  ++iview(i,4);
  BOOST_CHECK_EQUAL( aview(i,4),       5 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,4), 4 );
  BOOST_CHECK_EQUAL( int(iview(i,4)) , 4 );

  int val = iview(i,4);
  BOOST_CHECK_EQUAL( val, 4 );

  val = iview(i,4)+1;
  BOOST_CHECK_EQUAL( val, 5 );
}

BOOST_AUTO_TEST_CASE( test_indexview_3d )
{
  array::ArrayT<int> array( 5, 7, 10 );

  array::ArrayView<int,3>       aview(array);
  IndexView<int,3>       iview(array);
  const IndexView<int,3> const_iview(array);

  int i=2;
  int k=4;

  aview(i,0,k) = 1 IN_FORTRAN;
  BOOST_CHECK_EQUAL( aview(i,0,k),       1 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,0,k), 0 );
  BOOST_CHECK_EQUAL( int(iview(i,0,k)),  0 );

  iview(i,2,k) = 2;
  BOOST_CHECK_EQUAL( aview(i,2,k),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,2,k), 2 );
  BOOST_CHECK_EQUAL( int(iview(i,2,k)),  2 );

  iview(i,3,k) = iview(i,2,k);
  BOOST_CHECK_EQUAL( aview(i,3,k),       3 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,3,k), 2 );
  BOOST_CHECK_EQUAL( int(iview(i,3,k)),  2 );

  iview(i,3,k) = iview(i,2,k) + 1;
  BOOST_CHECK_EQUAL( aview(i,3,k),       4 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,3,k), 3 );
  BOOST_CHECK_EQUAL( int(iview(i,3,k)) , 3 );

  iview(i,4,k) = iview(i,3,k);
  ++iview(i,4,k);
  BOOST_CHECK_EQUAL( aview(i,4,k),       5 IN_FORTRAN );
  BOOST_CHECK_EQUAL( const_iview(i,4,k), 4 );
  BOOST_CHECK_EQUAL( int(iview(i,4,k)) , 4 );

  int val = iview(i,4,k);
  BOOST_CHECK_EQUAL( val, 4 );

  val = iview(i,4,k)+1;
  BOOST_CHECK_EQUAL( val, 5 );
}
*/

} // namespace test
} // namespace atlas
