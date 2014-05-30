/*
 * (C) Copyright 1996-2014 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestIndexView
#include <boost/test/included/unit_test.hpp>
#include "atlas/mpl/MPL.hpp"
#include "atlas/mesh/Array.hpp"
#include "atlas/mesh/ArrayView.hpp"
#include "atlas/mesh/IndexView.hpp"

using namespace atlas;

#ifdef HAVE_FORTRAN_NUMBERING
#define IN_FORTRAN
#else
#define IN_FORTRAN -1
#endif


BOOST_AUTO_TEST_CASE( test_indexview_1d )
{
  Array<int> array( 10 );

  ArrayView<int,1>       aview(array);
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
}

BOOST_AUTO_TEST_CASE( test_indexview_2d )
{
  Array<int> array( 5, 10 );

  ArrayView<int,2>       aview(array);
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
  Array<int> array( 5, 7, 10 );

  ArrayView<int,3>       aview(array);
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


