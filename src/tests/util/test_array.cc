/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestIndexView
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/Array.h"

using namespace atlas::array;

namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( test_array )
{
   auto ds = Array::create_storage<double>(4);
   auto hv = Array::make_host_view<double, 1>(ds);
   hv(3) = 4.5;

   BOOST_CHECK_EQUAL( hv(3) , 4.5 );
   
}

BOOST_AUTO_TEST_CASE( test_array_shape )
{
    ArrayShape as{2,3};
   auto ds = Array::create<double>(as);
   auto hv = Array::make_host_view<double, 2>(ds);
   hv(1,1) = 4.5;

   BOOST_CHECK_EQUAL( hv(1,1) , 4.5 );

}

}
}
