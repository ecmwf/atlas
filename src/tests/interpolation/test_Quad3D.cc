/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#define BOOST_TEST_MODULE atlas_test_geometry
#include "ecbuild/boost_test_framework.h"

#include "eckit/geometry/Point3.h"
#include "atlas/interpolation/Ray.h"
#include "atlas/interpolation/Quad3D.h"

using eckit::geometry::Point3;
using atlas::interpolation::Quad3D;
using atlas::interpolation::Intersect;
using atlas::interpolation::Ray;

//----------------------------------------------------------------------------------------------------------------------

const double relative_error = 0.0001;

BOOST_AUTO_TEST_CASE( test_quad_area )
{
    Point3 v0(0.,0.,0.);
    Point3 v1(1.,0.,0.);
    Point3 v2(1.,1.,0.);
    Point3 v3(0.,1.,0.);

    Quad3D quad1(v0.data(),v1.data(),v2.data(),v3.data());

    BOOST_CHECK( quad1.validate() );

    double area = quad1.area();

    BOOST_CHECK_CLOSE( area, 1.0, relative_error );

    Point3 c0(-2.,-2.,3.); // 4
    Point3 c1( 3.,-2.,3.); // 6
    Point3 c2( 3.,0.5,3.); // 1.5
    Point3 c3(-2.,0.5,3.); // 1

    Quad3D quad2(c0.data(),c1.data(),c2.data(),c3.data());

    BOOST_CHECK( quad2.validate() );

    area = quad2.area();

    BOOST_CHECK_CLOSE( area, 12.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_refquad )
{
  Point3 v0(0.,0.,0.);
  Point3 v1(1.,0.,0.);
  Point3 v2(1.,1.,0.);
  Point3 v3(0.,1.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.25,0.25,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.25, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.25, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_doublequad )
{
  Point3 v0(0.,0.,0.);
  Point3 v1(2.,0.,0.);
  Point3 v2(2.,2.,0.);
  Point3 v3(0.,2.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.5,0.5,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.25, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.25, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_rotatedquad )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.,0.,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.5, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_slopequad )
{
  Point3 v0( 2., 0.,2.);
  Point3 v1( 2., 0.,0.);
  Point3 v2( 0., 2.,0.);
  Point3 v3( 0., 2.,2.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(2.,2.,1.);
  Point3 dir (-1.,-1.,0.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.5, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_nointersect )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(2.,2.,1.);
  Point3 dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);
  BOOST_CHECK( ! isect );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_nointersect_aimoff )
{
  Point3 v0( 0.,-1.,0.);
  Point3 v1( 1., 0.,0.);
  Point3 v2( 0., 1.,0.);
  Point3 v3(-1., 0.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  Point3 orig(0.,0.,1.);
  Point3 dir (0.,1.,0.); // aim off

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);
  BOOST_CHECK( ! isect );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_corners )
{
  Point3 v0( 0.0,-2.0, 0.);
  Point3 v1( 2.5, 0.0, 0.);
  Point3 v2( 0.0, 3.5, 0.);
  Point3 v3(-1.5, 0.0, 0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  std::vector<Point3> corners;
  corners.push_back( Point3( 0.0,-2.0, 1.) );
  corners.push_back( Point3( 2.5, 0.0, 1.) );
  corners.push_back( Point3( 0.0, 3.5, 1.) );
  corners.push_back( Point3(-1.5, 0.0, 1.) );

  std::vector< std::pair<double,double> > uvs;
  uvs.push_back( std::make_pair(0.,0.) );
  uvs.push_back( std::make_pair(1.,0.) );
  uvs.push_back( std::make_pair(1.,1.) );
  uvs.push_back( std::make_pair(0.,1.) );

  for( size_t i = 0; i < 4; ++i )
  {
      Point3 orig = corners[i];
      Point3 dir (0.,0.,-1.);

      Ray ray(orig.data(),dir.data());

      Intersect isect = quad.intersects(ray);

      BOOST_CHECK( isect );
      BOOST_CHECK_CLOSE( isect.u, uvs[i].first , relative_error );
      BOOST_CHECK_CLOSE( isect.v, uvs[i].second, relative_error );
  }
}

