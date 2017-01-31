/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestFlags

#include <iostream>

#include "ecbuild/boost_test_framework.h"

#include "atlas/internals/Bitflags.h"
#include "tests/AtlasFixture.h"

using atlas::internals::Bitflags;

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_Flags )
{
  int b1 = (1<<0);
  int b2 = (1<<1);
  int b3 = (1<<2);
  int b4 = (1<<3);

	int bits = b1 | b2;
	std::cout << Bitflags::bitstr(bits) << std::endl;
	BOOST_CHECK_EQUAL( bits , 3);

	BOOST_CHECK_EQUAL( Bitflags::check(bits, b1 ) , true );
	BOOST_CHECK_EQUAL( Bitflags::check(bits, b2 ) , true );
	BOOST_CHECK_EQUAL( Bitflags::check(bits, b3 ) , false );
	BOOST_CHECK_EQUAL( Bitflags::check_all(bits, b1|b2 ) , true );
	BOOST_CHECK_EQUAL( Bitflags::check_all(bits, b1|b3 ) , false );
	BOOST_CHECK_EQUAL( Bitflags::check_any(bits, b1|b3 ) , true );
	BOOST_CHECK_EQUAL( Bitflags::check_any(bits, b3|b4 ) , false );
}

} // namespace test
} // namespace atlas
