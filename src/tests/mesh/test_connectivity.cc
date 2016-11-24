/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestConnectivity
#include "ecbuild/boost_test_framework.h"
#include "atlas/mesh/Connectivity.h"

using namespace atlas::mesh;

namespace atlas {
namespace mesh {

BOOST_AUTO_TEST_CASE( test_irregular_connectivity )
{
    IrregularConnectivity conn("mesh");
    BOOST_CHECK_EQUAL(conn.rows(),0);

    constexpr idx_t vals[4] = {2,3,5,6};
    conn.add(1, 4, vals, false);

    BOOST_CHECK_EQUAL(conn.rows(),1);
    BOOST_CHECK_EQUAL(conn.cols(0),4);

    BOOST_CHECK_EQUAL(conn.row(0)(0),2);
    BOOST_CHECK_EQUAL(conn.row(0)(1),3);
    BOOST_CHECK_EQUAL(conn.row(0)(2),5);
    BOOST_CHECK_EQUAL(conn.row(0)(3),6);

}


}
}
