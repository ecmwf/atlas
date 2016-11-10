/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestArray
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"

using namespace atlas::array;

namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( test_array )
{
   Array *ds = Array::create<double>(4ul);
   auto hv = make_gt_host_view<double, 1>(*ds);
   hv(3) = 4.5;

   ArrayView<double, 1> atlas_hv = make_host_view<double, 1>(*ds);

   BOOST_CHECK_EQUAL( hv(3) , 4.5 );
   BOOST_CHECK_EQUAL( atlas_hv(3) , 4.5 );

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_make_view )
{
   Array *ds = Array::create<double>(4ul);
   auto hv = make_gt_host_view<double, 1>(*ds);
   hv(3) = 4.5;

   ArrayView<double, 1> atlas_hv = make_view<double, 1>(*ds);

   BOOST_CHECK_EQUAL( hv(3) , 4.5 );
   BOOST_CHECK_EQUAL( atlas_hv(3) , 4.5 );
   
   delete ds;
}

BOOST_AUTO_TEST_CASE( test_localview )
{
   Array *ds = Array::create<double>(8ul,4ul,2ul);
   auto hv = make_view<double, 3>(*ds);
   
   // Initialize fields
   for( size_t i=0; i<ds->shape(0); ++i ) {
     for( size_t j=0; j<ds->shape(1); ++j ) {
       for( size_t k=0; k<ds->shape(2); ++k ) {
         hv(i,j,k) = (i*100) + (j*10) + (k);
       }
     }
   }

   // Check values
   for( size_t i=0; i<ds->shape(0); ++i ) {
     LocalView<double,2> lv = hv.at(0);
     for( size_t j=0; j<lv.shape(0); ++j ) {
       for( size_t k=0; k<lv.shape(1); ++k ) {
         BOOST_CHECK_EQUAL( lv(j,k), (i*100) + (j*10) + (k) );
       }
     }
   }
   
   delete ds;
}

BOOST_AUTO_TEST_CASE( test_array_shape )
{
   ArrayShape as{2,3};
   Array *ds = Array::create<double>(as);
   auto hv = make_gt_host_view<double, 2>(*ds);
   ArrayView<double, 2> atlas_hv = make_host_view<double, 2>(*ds);

   hv(1,1) = 4.5;

   BOOST_CHECK_EQUAL( hv(1,1) , 4.5 );
   BOOST_CHECK_EQUAL( atlas_hv(1,1) , 4.5 );

   BOOST_CHECK_EQUAL( ds->size() , 6 );
   BOOST_CHECK_EQUAL( ds->rank() , 2 );
   BOOST_CHECK_EQUAL( ds->stride(0) , 3 );
   BOOST_CHECK_EQUAL( ds->stride(1) , 1 );
   BOOST_CHECK_EQUAL( ds->contiguous() , true );

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_spec )
{
   Array* ds = Array::create<double>(4,5,6);
   BOOST_CHECK_EQUAL( ds->spec().rank(), 3);
   BOOST_CHECK_EQUAL( ds->spec().size(), 4*5*6);
   BOOST_CHECK_EQUAL( ds->spec().shape()[0], 4);
   BOOST_CHECK_EQUAL( ds->spec().shape()[1], 5);
   BOOST_CHECK_EQUAL( ds->spec().shape()[2], 6);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[0], 6);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[1], 5);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[2], 4);

   BOOST_CHECK_EQUAL( ds->spec().strides()[0], 6*5);
   BOOST_CHECK_EQUAL( ds->spec().strides()[1], 6);
   BOOST_CHECK_EQUAL( ds->spec().strides()[2], 1);

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_spec_layout )
{
   Array* ds = Array::create_with_layout<double, gridtools::layout_map<0,1,2> >(4,5,6);
   BOOST_CHECK_EQUAL( ds->spec().rank(), 3);
   BOOST_CHECK_EQUAL( ds->spec().size(), 4*5*6);
   BOOST_CHECK_EQUAL( ds->spec().shape()[0], 4);
   BOOST_CHECK_EQUAL( ds->spec().shape()[1], 5);
   BOOST_CHECK_EQUAL( ds->spec().shape()[2], 6);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[0], 6);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[1], 5);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[2], 4);
   BOOST_CHECK_EQUAL( ds->spec().strides()[0], 6*5);
   BOOST_CHECK_EQUAL( ds->spec().strides()[1], 6);
   BOOST_CHECK_EQUAL( ds->spec().strides()[2], 1);

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_spec_layout_rev )
{
   Array* ds = Array::create_with_layout<double, gridtools::layout_map<2,1,0> >(4,5,6);
   BOOST_CHECK_EQUAL( ds->spec().rank(), 3);
   BOOST_CHECK_EQUAL( ds->spec().size(), 4*5*6);
   BOOST_CHECK_EQUAL( ds->spec().shape()[0], 4);
   BOOST_CHECK_EQUAL( ds->spec().shape()[1], 5);
   BOOST_CHECK_EQUAL( ds->spec().shape()[2], 6);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[0], 4);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[1], 5);
   BOOST_CHECK_EQUAL( ds->spec().shapef()[2], 6);
   BOOST_CHECK_EQUAL( ds->spec().strides()[0], 1);
   BOOST_CHECK_EQUAL( ds->spec().strides()[1], 4);
   BOOST_CHECK_EQUAL( ds->spec().strides()[2], 4*5);

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_resize_throw )
{
   Array* ds = Array::create<double>(32,5,33);

   BOOST_CHECK_NO_THROW(ds->resize(32,5,33));
   BOOST_CHECK_THROW(ds->resize(32,4,33), eckit::BadParameter);
   BOOST_CHECK_THROW(ds->resize(32,5,32), eckit::BadParameter);
   BOOST_CHECK_THROW(ds->resize(32,5,33,4), eckit::BadParameter);

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_resize )
{
  Array* ds = Array::create<double>(7, 5, 8);
  {
    ArrayView<double, 3> hv = make_host_view<double, 3>(*ds);
    hv(3, 3, 3) = 4.5;
    hv(6, 4, 7) = 7.5;
  }
  ds->resize(32, 5, 33);
  ArrayView<double, 3> hv = make_host_view<double, 3>(*ds);

  BOOST_CHECK_EQUAL(ds->spec().shape()[0], 32);
  BOOST_CHECK_EQUAL(ds->spec().shape()[1], 5);
  BOOST_CHECK_EQUAL(ds->spec().shape()[2], 33);

  BOOST_CHECK_EQUAL(ds->spec().rank(), 3);
  BOOST_CHECK_EQUAL(ds->spec().size(), 32 * 5 * 33);

  BOOST_CHECK_EQUAL(hv(3, 3, 3), 4.5);
  BOOST_CHECK_EQUAL(hv(6, 4, 7), 7.5);

  delete ds;
}

BOOST_AUTO_TEST_CASE( test_resize_shape )
{
   Array* ds = Array::create<double>(7,5,8);
   {
     ArrayView<double, 3> hv = make_host_view<double, 3>(*ds);
     hv(3, 3, 3) = 4.5;
     hv(6, 4, 7) = 7.5;
   }
   ds->resize(ArrayShape{32,5,33});

   ArrayView<double, 3> hv = make_host_view<double, 3>(*ds);
   BOOST_CHECK_EQUAL(ds->spec().shape()[0], 32);
   BOOST_CHECK_EQUAL(ds->spec().shape()[1], 5);
   BOOST_CHECK_EQUAL(ds->spec().shape()[2], 33);

   BOOST_CHECK_EQUAL( ds->spec().rank(), 3);
   BOOST_CHECK_EQUAL( ds->spec().size(), 32*5*33);

   BOOST_CHECK_EQUAL( hv(3,3,3), 4.5);
   BOOST_CHECK_EQUAL( hv(6,4,7), 7.5);

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_insert )
{
   Array* ds = Array::create<double>(7,5,8);

   ArrayView<double, 3> hv = make_host_view<double, 3>(*ds);
   hv(1,3,3) = 1.5;
   hv(2,3,3) = 2.5;
   hv(3,3,3) = 3.5;
   hv(6,4,7) = 6.5;

   ds->insert(2, 2);

   BOOST_CHECK_EQUAL(ds->spec().shape()[0], 9);
   BOOST_CHECK_EQUAL(ds->spec().shape()[1], 5);
   BOOST_CHECK_EQUAL(ds->spec().shape()[2], 8);

   BOOST_CHECK_EQUAL( ds->spec().rank(), 3);
   BOOST_CHECK_EQUAL( ds->spec().size(), 9*5*8);

   ArrayView<double, 3> hv2 = make_host_view<double, 3>(*ds);

   BOOST_CHECK_EQUAL( hv(1,3,3), 1.5);
   BOOST_CHECK_EQUAL( hv(3,3,3), 3.5);

   BOOST_CHECK_EQUAL( hv2(1,3,3), 1.5);
   BOOST_CHECK_EQUAL( hv2(2,3,3), 2.5);
   BOOST_CHECK_EQUAL( hv2(3,3,3), 3.5);
   BOOST_CHECK_EQUAL( hv2(6,4,7), 6.5);

   delete ds;
}

BOOST_AUTO_TEST_CASE( test_wrap_storage )
{
   Array* ds = Array::create<double>(4,5,6);

   ArrayView<double, 3> hv = make_host_view<double, 3>(*ds);

   hv(2,3,3) = 2.5;
   Array* ds_ext = Array::wrap<double>(hv.data(), ds->spec());

   //TODO seg fault
   ArrayView<double, 3> hv_ext = make_host_view<double, 3>(*ds_ext);

//   BOOST_CHECK_EQUAL( hv_ext(2,3,3), 2.5);

//   delete ds;
//   delete ds_ext;
}

}
}
