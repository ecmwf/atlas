/*
 * (C) Copyright 1996-2012 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <string>
#include <iostream>

#define BOOST_TEST_MODULE TestRotatedLatLon
#define BOOST_UNIT_TEST_FRAMEWORK_HEADER_ONLY
#include "ecbuild/boost_test_framework.h"

#include "eckit/filesystem/PathName.h"
#include "eckit/filesystem/LocalPathName.h"

#include "eckit/grib/GribHandle.h"
#include "eckit/grib/GribAccessor.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/Grib.h"
#include "atlas/grid/GridSpec.h"
#include "atlas/grid/RotatedLatLon.h"


using namespace std;
using namespace eckit;
using namespace eckit::grib;
using namespace atlas;
using namespace atlas::grid;


BOOST_AUTO_TEST_SUITE( TestRotatedLatLon )


BOOST_AUTO_TEST_CASE( test_rotated_lat_lon )
{
   cout << "\nGrid:: ...test_rotated_lat_lon\n";

   Grid::Point south_pole(37.5,177.5);
   double polerot = 10.0 ;

   Grid::Point point(51.0,-3.0);

   Rotgrid mapping(south_pole, polerot);
   Grid::Point tr_point = mapping.transform(point);
   Grid::Point point2 = mapping.transform(tr_point, true);

   std::cout << " Location of pole of rotated grid as seen in non-rotated grid: " << south_pole << "\n";
   std::cout << " Additional axial rotation about pole of rotated grid: polerot= " << polerot << "\n";
   std::cout << " Location of chosen point in non-rotated grid           : " << point << "\n";
   std::cout << " Location of chosen point as seen in rotated grid       : " << tr_point << "\n";
   std::cout << " Location of chosen point put back into non-rotated grid: " << point2 << "\n";
   BOOST_CHECK_CLOSE(point2.lat(),point.lat(),Grid::degrees_eps());
   BOOST_CHECK_CLOSE(point2.lon(),point.lon(),Grid::degrees_eps());

   {
       Grid::Point tr_point = mapping.transform(point);
       Grid::Point point2 = mapping.transform(tr_point);
       std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " rotate        = " << tr_point << " unrotated " << point2 << "\n";

       Grid::Point rotated = mapping.magics_rotate(point);
       Grid::Point unrotated = mapping.magics_rotate(rotated);
       std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " Magic rotated = " <<  rotated << " unrotated " << unrotated << "\n";

       Grid::Point rotated1 = mapping.rotate(point);
       Grid::Point unrotated2 = mapping.unrotate(rotated1);
       std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " My    rotated = " <<  rotated1 << " unrotated " << unrotated2 << "\n";
    }
}


BOOST_AUTO_TEST_CASE( test_rotated_lat_lon_2 )
{
   cout << "\nGrid:: ...test_rotated_lat_lon_2\n";

   // This was taken form the web :http://www.mathworks.se/matlabcentral/fileexchange/43435-rotated-grid-transform
   // Used to compare with our transform
   //
   //This functions transforms a set of coordinates in regular lon/lat degrees,
   // grid_in = [lon, lat], to a set of coordinates in rotated lon/lat degrees, grid_out = [lon', lat'], and vice versa:
   //
   //[grid_out] = rotated_grid_transform(grid_in, option, SP_coor)
   //
   //where option is the 'direction' of the transform (1: regular -> rotated and 2: rotated -> regular)
   //and SP_coor are the coordinates of the South Pole in the rotated grid [SP_lon, SP_lat]
   //
   //Example:
   //SP_coor = [18 -39.3];
   //grid_in = [[12; 12; 12],[55; 54; 53]];
   //[grid_out] = rotated_grid_transform(grid_in, 1, SP_coor)
   //
   //grid_out =
   //
   //   -3.4476 4.4397
   //   -3.5289 3.4430
   //   -3.6100 2.4463
   //
   //grid_in = grid_out;
   //[grid_out] = rotated_grid_transform(grid_in, 2, SP_coor)
   //
   //grid_out =
   //
   //   12.0000 55.0000
   //   12.0000 54.0000
   //   12.0000 53.0000

   Grid::Point south_pole(18,-39.3);
   double polerot = 0.0 ;
   Rotgrid mapping(south_pole, polerot);
   {
      Grid::Point point(12,55);
      Grid::Point tr_point = mapping.transform(point);
      Grid::Point point2 = mapping.transform(tr_point);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " rotate        = " << tr_point << " unrotated " << point2 << "\n";

      Grid::Point rotated = mapping.magics_rotate(point);
      Grid::Point unrotated = mapping.magics_rotate(rotated);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " Magic rotated = " <<  rotated << " unrotated " << unrotated << "\n";

      Grid::Point rotated1 = mapping.rotate(point);
      Grid::Point unrotated2 = mapping.unrotate(rotated1);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " My    rotated = " <<  rotated1 << " unrotated " << unrotated2 << "\n";
      BOOST_CHECK_CLOSE(point.lat(),unrotated2.lat(),Grid::degrees_eps());
      BOOST_CHECK_CLOSE(point.lon(),unrotated2.lon(),Grid::degrees_eps());
   }

   {
      Grid::Point point(12,54);
      Grid::Point tr_point = mapping.transform(point);
      Grid::Point point2 = mapping.transform(tr_point);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " rotate        = " << tr_point << " unrotated " << point2 << "\n";

      Grid::Point rotated = mapping.magics_rotate(point);
      Grid::Point unrotated = mapping.magics_rotate(rotated);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " Magic rotated = " <<  rotated << " unrotated " << unrotated << "\n";

      Grid::Point rotated1 = mapping.rotate(point);
      Grid::Point unrotated2 = mapping.unrotate(rotated1);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " My    rotated = " <<  rotated1 << " unrotated " << unrotated2 << "\n";
      BOOST_CHECK_CLOSE(point.lat(),unrotated2.lat(),Grid::degrees_eps());
      BOOST_CHECK_CLOSE(point.lon(),unrotated2.lon(),Grid::degrees_eps());
   }

   {
      Grid::Point point(12,53);
      Grid::Point tr_point = mapping.transform(point);
      Grid::Point point2 = mapping.transform(tr_point);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " rotate        = " << tr_point << " unrotated " << point2 << "\n";

      Grid::Point rotated = mapping.magics_rotate(point);
      Grid::Point unrotated = mapping.magics_rotate(rotated);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " Magic rotated = " <<  rotated << " unrotated " << unrotated << "\n";

      Grid::Point rotated1 = mapping.rotate(point);
      Grid::Point unrotated2 = mapping.unrotate(rotated1);
      std::cout <<  " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " My    rotated = " <<  rotated1 << " unrotated " << unrotated2 << "\n";
      BOOST_CHECK_CLOSE(point.lat(),unrotated2.lat(),Grid::degrees_eps());
      BOOST_CHECK_CLOSE(point.lon(),unrotated2.lon(),Grid::degrees_eps());
   }
}

BOOST_AUTO_TEST_CASE( test_south_pole_at_minus_90_0 )
{
   cout << "\nGrid:: ...test_south_pole_at_minus_90_0\n";

   // Rotation of -90,0 for south pole, should mean no change, since -90,0 is at the south pole
   // Should end up with identity matrix and hence no change, likewise for un-rotate
   Grid::Point south_pole(-90,0);
   double polerot = 0.0 ;
   Rotgrid mapping(south_pole, polerot);
   {
      Grid::Point point(12,55);
      Grid::Point tr_point = mapping.transform(point);
      Grid::Point point2 = mapping.transform(tr_point);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " rotate        = " << tr_point << " unrotated " << point2 << "\n";

      Grid::Point rotated = mapping.magics_rotate(point);
      Grid::Point unrotated = mapping.magics_rotate(rotated);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " Magic rotated = " <<  rotated << " unrotated " << unrotated << "\n";

      Grid::Point rotated1 = mapping.rotate(point);
      Grid::Point unrotated2 = mapping.unrotate(rotated1);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " My    rotated = " <<  rotated1 << " unrotated " << unrotated2 << "\n";
      BOOST_CHECK_CLOSE(point.lat(),rotated1.lat(),Grid::degrees_eps());
      BOOST_CHECK_CLOSE(point.lon(),rotated1.lon(),Grid::degrees_eps());
      BOOST_CHECK_CLOSE(point.lat(),unrotated2.lat(),Grid::degrees_eps());
      BOOST_CHECK_CLOSE(point.lon(),unrotated2.lon(),Grid::degrees_eps());
   }
}

BOOST_AUTO_TEST_CASE( test_south_pole_at_minus_0_0 )
{
   cout << "\nGrid:: ...test_south_pole_at_minus_0_0\n";

   // Move south pole to the equator
   // Henve a poit at the south pole ,  shold move to the equator
   Grid::Point south_pole(0,0); double polerot = 0.0 ;
   Rotgrid mapping(south_pole, polerot);
   {
      Grid::Point point(-90,0);
      Grid::Point tr_point = mapping.transform(point);
      Grid::Point point2 = mapping.transform(tr_point);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " rotate        = " << tr_point << " unrotated " << point2 << "\n";

      Grid::Point rotated = mapping.magics_rotate(point);
      Grid::Point unrotated = mapping.magics_rotate(rotated);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " Magic rotated = " <<  rotated << " unrotated " << unrotated << "\n";

      Grid::Point rotated1 = mapping.rotate(point);
      Grid::Point unrotated2 = mapping.unrotate(rotated1);
      std::cout << " sp" << south_pole << " sp_rot(" << polerot << ") " << point << " My    rotated = " <<  rotated1 << " unrotated " << unrotated2 << "\n";
   }
}


BOOST_AUTO_TEST_SUITE_END()
