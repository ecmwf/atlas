/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#define BOOST_TEST_MODULE TestPolygon
#include "ecbuild/boost_test_framework.h"
#include "tests/AtlasFixture.h"

#include <algorithm>
#include <cmath>
#include "atlas/mesh/detail/Polygon.h"

using atlas::mesh::detail::Polygon;

namespace atlas {
namespace test {

BOOST_GLOBAL_FIXTURE( AtlasFixture );

BOOST_AUTO_TEST_SUITE( test_polygon )

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( test_polygon_coordinates )

BOOST_AUTO_TEST_CASE( test_polygon_something )
{
    struct : Polygon {
        void setCoordinates() const {
            includesNorthPole_ = includesSouthPole_ = false;

            using p = PointLonLat;
            coordinates_ = {
                p(122.143, 35.9951),
                p(120, 30.4576),
                p(118.125, 24.9199),
                p(116.471, 19.3822),
                p(121.765, 19.3822),
                p(120, 13.8445),
                p(118.421, 8.3067),
                p(117, 2.7689),
                p(112.5, 2.7689),
                p(112.5, -2.7689),
                p(117, -2.7689),
                p(113.684, -8.3067),
                p(118.421, -8.3067),
                p(115, -13.8445),
                p(120, -13.8445),
                p(116.471, -19.3822),
                p(112.5, -24.9199),
                p(118.125, -24.9199),
                p(114, -30.4576),
                p(120, -30.4576),
                p(115.714, -35.9951),
                p(122.143, -35.9951),
                p(128.571, -35.9951),
                p(135, -35.9951),
                p(141.429, -35.9951),
                p(147.857, -35.9951),
                p(154.286, -35.9951),
                p(160.714, -35.9951),
                p(167.143, -35.9951),
                p(173.571, -35.9951),
                p(180, -35.9951),
                p(186.429, -35.9951),
                p(192.857, -35.9951),
                p(199.286, -35.9951),
                p(205.714, -35.9951),
                p(212.143, -35.9951),
                p(218.571, -35.9951),
                p(225, -35.9951),
                p(231.429, -35.9951),
                p(237.857, -35.9951),
                p(240, -30.4576),
                p(234, -30.4576),
                p(236.25, -24.9199),
                p(238.235, -19.3822),
                p(232.941, -19.3822),
                p(235, -13.8445),
                p(236.842, -8.3067),
                p(238.5, -2.7689),
                p(234, -2.7689),
                p(234, 2.7689),
                p(238.5, 2.7689),
                p(241.579, 8.3067),
                p(236.842, 8.3067),
                p(240, 13.8445),
                p(235, 13.8445),
                p(238.235, 19.3822),
                p(241.875, 24.9199),
                p(236.25, 24.9199),
                p(240, 30.4576),
                p(244.286, 35.9951),
                p(237.857, 35.9951),
                p(231.429, 35.9951),
                p(225, 35.9951),
                p(218.571, 35.9951),
                p(212.143, 35.9951),
                p(205.714, 35.9951),
                p(199.286, 35.9951),
                p(192.857, 35.9951),
                p(186.429, 35.9951),
                p(180, 35.9951),
                p(173.571, 35.9951),
                p(167.143, 35.9951),
                p(160.714, 35.9951),
                p(154.286, 35.9951),
                p(147.857, 35.9951),
                p(141.429, 35.9951),
                p(135, 35.9951),
                p(128.571, 35.9951),
                p(122.143, 35.9951)
            };

            coordinatesMin_ = coordinatesMax_ = coordinates_.front();
            for (const PointLonLat& pt : coordinates_) {
                coordinatesMin_ = PointLonLat::componentsMin(coordinatesMin_, pt);
                coordinatesMax_ = PointLonLat::componentsMax(coordinatesMax_, pt);
            }
        }

    } poly;

    std::vector<idx_t> indices(79);
    std::iota(indices.begin(), indices.end(), 0);
    for (auto i : indices) {
        poly.push_back(i);
    }
    poly.setCoordinates();

    PointLonLat P[] = {
        PointLonLat(-118.8, 32.0919),
        PointLonLat(-118.929, 23.7202),
        PointLonLat(118.8, -32.0919)
    };

//    BOOST_CHECK_EQUAL(poly.containsPointInLonLatGeometry(P[2]), true);
    BOOST_CHECK_EQUAL(poly.containsPointInSphericalGeometry(P[2]), false);
}

BOOST_AUTO_TEST_SUITE_END()

// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace atlas
