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
#include "atlas/array/Vector.h"

using namespace atlas::array;

namespace atlas {
namespace test {

BOOST_AUTO_TEST_CASE(test_vector) {
    Vector<int> vec(3);

    vec[0] = 3;
    vec[1] = -3;
    vec[2] = 1;

    BOOST_CHECK_EQUAL( vec.size(),3);
    BOOST_CHECK_EQUAL( vec[0],3);
    BOOST_CHECK_EQUAL( vec[1],-3);
    BOOST_CHECK_EQUAL( vec[2],1);

    vec.resize(5);
    vec[3] = 5;
    vec[4] = 6;
    BOOST_CHECK_EQUAL( vec[0],3);
    BOOST_CHECK_EQUAL( vec[1],-3);
    BOOST_CHECK_EQUAL( vec[2],1);
    BOOST_CHECK_EQUAL( vec[3],5);
    BOOST_CHECK_EQUAL( vec[4],6);

}

}
}
