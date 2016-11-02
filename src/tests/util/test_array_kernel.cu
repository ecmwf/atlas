/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestArrayKernel
#include "ecbuild/boost_test_framework.h"
#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"

using namespace atlas::array;

namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE( test_array )
{
   auto ds = Array::create<double>(4ul);
   auto hv = make_gt_host_view<double, 1>(ds);
   hv(3) = 4.5;

   ArrayView<double, 1> atlas_hv = make_host_view<double, 1>(ds);

   BOOST_CHECK_EQUAL( hv(3) , 4.5 );
   BOOST_CHECK_EQUAL( atlas_hv(3) , 4.5 );
}

}
}
