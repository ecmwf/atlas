/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <algorithm>
#include <iomanip>
#include <sstream>

#define BOOST_TEST_MODULE TestDomain
#include "ecbuild/boost_test_framework.h"


#include "atlas/domain.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "tests/AtlasFixture.h"

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_CASE( test_domain_rectangular )
{
  Domain domain = RectangularDomain( {0,180}, {-25,25} );
  BOOST_CHECK( not domain.global() );
  BOOST_CHECK_EQUAL( domain.type(), std::string("rectangular") );

  util::Config domain_cfg = domain.spec();
  Domain from_cfg(domain_cfg);
  Log::info() << from_cfg.spec() << std::endl;
  BOOST_CHECK_EQUAL( from_cfg.type(), std::string("rectangular") );
}

BOOST_AUTO_TEST_CASE( test_domain_rectangular_tolerance )
{
    using grid::StructuredGrid;
    using grid::LinearSpacing;

    Domain domain = RectangularDomain( {-27, 45}, {33, 73} );
    RectangularDomain rd = domain;

    StructuredGrid::XSpace xspace( LinearSpacing( rd.xmin(), rd.xmax(), 721, true ));
    StructuredGrid::YSpace yspace( LinearSpacing( rd.ymin(), rd.ymax(), 401 ));

    StructuredGrid grid(xspace, yspace, StructuredGrid::Projection(), rd);

    BOOST_CHECK_EQUAL(grid.nx(0), 721);
    BOOST_CHECK_EQUAL(grid.ny(), 401);
}

BOOST_AUTO_TEST_CASE( test_domain_zonal_from_rectangular )
{
  Domain domain = RectangularDomain( {0,360}, {-25,25} );
  BOOST_CHECK( not domain.global() );
  BOOST_CHECK_EQUAL( domain.type(), std::string("zonal_band") );

  util::Config domain_cfg = domain.spec();
  Domain from_cfg(domain_cfg);
  Log::info() << from_cfg.spec() << std::endl;
  BOOST_CHECK_EQUAL( from_cfg.type(), std::string("zonal_band") );
}

BOOST_AUTO_TEST_CASE( test_domain_global_from_rectangular )
{
  Domain domain = RectangularDomain( {-180,180}, {-90,90} );
  BOOST_CHECK( domain.global() );
  BOOST_CHECK_EQUAL( domain.type(), std::string("global") );

  util::Config domain_cfg = domain.spec();
  Domain from_cfg(domain_cfg);
  Log::info() << from_cfg.spec() << std::endl;
  BOOST_CHECK_EQUAL( from_cfg.type(), std::string("global") );

  RectangularDomain rd = domain;
  BOOST_CHECK( rd == true );
  BOOST_CHECK_EQUAL( rd.xmin(), -180.);
  BOOST_CHECK_EQUAL( rd.xmax(),  180.);
  BOOST_CHECK_EQUAL( rd.ymin(), -90. );
  BOOST_CHECK_EQUAL( rd.ymax(),  90. );

  ZonalBandDomain zd = domain;
  BOOST_CHECK( zd == true );
  BOOST_CHECK_EQUAL( zd.xmin(), -180.);
  BOOST_CHECK_EQUAL( zd.xmax(),  180.);
  BOOST_CHECK_EQUAL( zd.ymin(), -90. );
  BOOST_CHECK_EQUAL( zd.ymax(),  90. );
}

BOOST_AUTO_TEST_CASE( test_domain_global_from_zonalband )
{
  Domain domain = ZonalBandDomain( {-45,45} );
  BOOST_CHECK( not domain.global() );
  BOOST_CHECK_EQUAL( domain.type(), std::string("zonal_band") );

  util::Config domain_cfg = domain.spec();
  Domain from_cfg(domain_cfg);
  Log::info() << from_cfg.spec() << std::endl;
  BOOST_CHECK_EQUAL( from_cfg.type(), std::string("zonal_band") );

  RectangularDomain rd = domain;
  BOOST_CHECK( rd == true );
  BOOST_CHECK_EQUAL( rd.xmin(), 0);
  BOOST_CHECK_EQUAL( rd.xmax(), 360.);
  BOOST_CHECK_EQUAL( rd.ymin(), -45. );
  BOOST_CHECK_EQUAL( rd.ymax(),  45. );

  ZonalBandDomain zd = domain;
  BOOST_CHECK( zd == true );
  BOOST_CHECK_EQUAL( zd.ymin(), -45. );
  BOOST_CHECK_EQUAL( zd.ymax(),  45. );
}


} // namespace test
} // namespace atlas
