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

#include "atlas/interpolation/method/Ray.h"
#include "atlas/interpolation/element/Quad3D.h"
#include "atlas/util/Point.h"

using atlas::interpolation::element::Quad3D;
using atlas::interpolation::method::Intersect;
using atlas::interpolation::method::Ray;
using atlas::PointXYZ;

//----------------------------------------------------------------------------------------------------------------------

const double relative_error = 0.0001;

BOOST_AUTO_TEST_CASE( test_quad_area )
{
    PointXYZ v0(0.,0.,0.);
    PointXYZ v1(1.,0.,0.);
    PointXYZ v2(1.,1.,0.);
    PointXYZ v3(0.,1.,0.);

    Quad3D quad1(v0.data(),v1.data(),v2.data(),v3.data());

    BOOST_CHECK( quad1.validate() );

    double area = quad1.area();

    BOOST_CHECK_CLOSE( area, 1.0, relative_error );

    PointXYZ c0(-2.,-2.,3.); // 4
    PointXYZ c1( 3.,-2.,3.); // 6
    PointXYZ c2( 3.,0.5,3.); // 1.5
    PointXYZ c3(-2.,0.5,3.); // 1

    Quad3D quad2(c0.data(),c1.data(),c2.data(),c3.data());

    BOOST_CHECK( quad2.validate() );

    area = quad2.area();

    BOOST_CHECK_CLOSE( area, 12.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_refquad )
{
  PointXYZ v0(0.,0.,0.);
  PointXYZ v1(1.,0.,0.);
  PointXYZ v2(1.,1.,0.);
  PointXYZ v3(0.,1.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  PointXYZ orig(0.25,0.25,1.);
  PointXYZ dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.25, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.25, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_doublequad )
{
  PointXYZ v0(0.,0.,0.);
  PointXYZ v1(2.,0.,0.);
  PointXYZ v2(2.,2.,0.);
  PointXYZ v3(0.,2.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  PointXYZ orig(0.5,0.5,1.);
  PointXYZ dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.25, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.25, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_rotatedquad )
{
  PointXYZ v0( 0.,-1.,0.);
  PointXYZ v1( 1., 0.,0.);
  PointXYZ v2( 0., 1.,0.);
  PointXYZ v3(-1., 0.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  PointXYZ orig(0.,0.,1.);
  PointXYZ dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.5, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_slopequad )
{
  PointXYZ v0( 2., 0.,2.);
  PointXYZ v1( 2., 0.,0.);
  PointXYZ v2( 0., 2.,0.);
  PointXYZ v3( 0., 2.,2.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  PointXYZ orig(2.,2.,1.);
  PointXYZ dir (-1.,-1.,0.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);

  BOOST_CHECK( isect );
  BOOST_CHECK_CLOSE( isect.u, 0.5, relative_error );
  BOOST_CHECK_CLOSE( isect.v, 0.5, relative_error );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_nointersect )
{
  PointXYZ v0( 0.,-1.,0.);
  PointXYZ v1( 1., 0.,0.);
  PointXYZ v2( 0., 1.,0.);
  PointXYZ v3(-1., 0.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  PointXYZ orig(2.,2.,1.);
  PointXYZ dir (0.,0.,-1.);

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);
  BOOST_CHECK( ! isect );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_nointersect_aimoff )
{
  PointXYZ v0( 0.,-1.,0.);
  PointXYZ v1( 1., 0.,0.);
  PointXYZ v2( 0., 1.,0.);
  PointXYZ v3(-1., 0.,0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  PointXYZ orig(0.,0.,1.);
  PointXYZ dir (0.,1.,0.); // aim off

  Ray ray(orig.data(),dir.data());

  Intersect isect = quad.intersects(ray);
  BOOST_CHECK( ! isect );
}

BOOST_AUTO_TEST_CASE( test_quadrilateral_intersection_corners )
{
  PointXYZ v0( 0.0,-2.0, 0.);
  PointXYZ v1( 2.5, 0.0, 0.);
  PointXYZ v2( 0.0, 3.5, 0.);
  PointXYZ v3(-1.5, 0.0, 0.);

  Quad3D quad(v0.data(),v1.data(),v2.data(),v3.data());

  BOOST_CHECK( quad.validate() );

  std::vector<PointXYZ> corners;
  corners.push_back( PointXYZ( 0.0,-2.0, 1.) );
  corners.push_back( PointXYZ( 2.5, 0.0, 1.) );
  corners.push_back( PointXYZ( 0.0, 3.5, 1.) );
  corners.push_back( PointXYZ(-1.5, 0.0, 1.) );

  std::vector< std::pair<double,double> > uvs;
  uvs.push_back( std::make_pair(0.,0.) );
  uvs.push_back( std::make_pair(1.,0.) );
  uvs.push_back( std::make_pair(1.,1.) );
  uvs.push_back( std::make_pair(0.,1.) );

  for( size_t i = 0; i < 4; ++i )
  {
      PointXYZ orig = corners[i];
      PointXYZ dir (0.,0.,-1.);

      Ray ray(orig.data(),dir.data());

      Intersect isect = quad.intersects(ray);

      BOOST_CHECK( isect );
      BOOST_CHECK_CLOSE( isect.u, uvs[i].first , relative_error );
      BOOST_CHECK_CLOSE( isect.v, uvs[i].second, relative_error );
  }
}

