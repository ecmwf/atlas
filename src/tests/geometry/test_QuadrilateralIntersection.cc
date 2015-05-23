/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#define BOOST_TEST_MODULE test_atlas_geometry
#include "ecbuild/boost_test_framework.h"

#include "eckit/geometry/Point3.h"
#include "atlas/geometry/QuadrilateralIntersection.h"

using eckit::geometry::Point3;
using eckit::geometry::points_equal;
using atlas::geometry::QuadrilateralIntersection;
using atlas::geometry::Intersect;
using atlas::geometry::Ray;

//----------------------------------------------------------------------------------------------------------------------

const double relative_error = 0.0001;

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_refquad )
{
  Point3 v0(0.,0.,0.);
  Point3 v1(1.,0.,0.);
  Point3 v2(1.,1.,0.);
  Point3 v3(0.,1.,0.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.25,0.25,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( quad.validateIntersection(ray) );

  BOOST_CHECK_CLOSE( isect.u, 0.25, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.25, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_doublequad )
{
  Point3 v0(0.,0.,0.);
  Point3 v1(2.,0.,0.);
  Point3 v2(2.,2.,0.);
  Point3 v3(0.,2.,0.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.5,0.5,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( quad.validateIntersection(ray) );

  BOOST_CHECK_CLOSE( isect.u, 0.25, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.25, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_rotatedquad )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.,0.,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( quad.validateIntersection(ray) );

  BOOST_CHECK_CLOSE( isect.u, 0.5, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_slopequad )
{
  Point3 v0( 2., 0.,2.);
  Point3 v1( 2., 0.,0.);
  Point3 v2( 0., 2.,0.);
  Point3 v3( 0., 2.,2.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(2.,2.,1.);
  Point3 dir (-1.,-1.,0.);

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( quad.validateIntersection(ray) );

  BOOST_CHECK_CLOSE( isect.u, 0.5, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_corner )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(-1.,0.,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( quad.validateIntersection(ray) );

  BOOST_CHECK_CLOSE( isect.u, 0.0, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 1.0, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_nointersect )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(2.,2.,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( ! isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( ! quad.validateIntersection(ray) );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_nointersect_aimoff )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  QuadrilateralIntersection quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.,0.,1.);
  Point3 dir (0.,1.,0.); // aim off

  Ray ray(orig.data(),dir.data());
  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( ! isect );

  BOOST_TEST_MESSAGE( isect );

  BOOST_CHECK( ! quad.validateIntersection(ray) );
}

