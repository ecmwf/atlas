/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <vector>

#include "atlas/runtime/Log.h"
#include "atlas/util/Config.h"
#include "atlas/util/Constants.h"
#include "atlas/util/Rotation.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::util::Config;
using atlas::util::Rotation;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

constexpr double eps = 1.e-5;
constexpr double d2r = atlas::util::Constants::degreesToRadians();
constexpr double r2d = atlas::util::Constants::radiansToDegrees();

bool equivalent(const PointLonLat& p1, const PointLonLat& p2) {
    using eckit::types::is_approximately_equal;
    auto f = [=](double lon) { return 10. + std::cos(lon * d2r); };

    return is_approximately_equal(p1.lat(), p2.lat(), eps) &&
           (is_approximately_equal(std::abs(p2.lat()), 90., eps) ||
            is_approximately_equal(f(p1.lon()), f(p2.lon()), eps));
}

#define EXPECT_EQUIVALENT(p1, p2) EXPECT(equivalent(p1, p2))

//-----------------------------------------------------------------------------

class MagicsRotation {
    // For reference, this what magics uses, it appears as if it originated from
    // fortran code
    // Strangely the definition of rotate and unrotate are switched.

public:
    MagicsRotation(const PointLonLat& south_pole): south_pole_(south_pole) {}

    PointLonLat rotate(const PointLonLat& point) const {
        return magics_unrotate(point);  /// Switch meaning !!!
    }

    PointLonLat unrotate(const PointLonLat& point) const {
        return magics_rotate(point);  /// Swich meaning !!!
    }

private:
    PointLonLat south_pole_;

    PointLonLat magics_rotate(const PointLonLat& point) const {
        double lat_y = point.lat();
        double lon_x = point.lon();

        double sin_south_pole_lat = std::sin(d2r * (south_pole_.lat() + 90.));
        double cos_south_pole_lat = std::cos(d2r * (south_pole_.lat() + 90.));

        double ZXMXC           = d2r * (lon_x - south_pole_.lon());
        double sin_lon_decr_sp = std::sin(ZXMXC);
        double cos_lon_decr_sp = std::cos(ZXMXC);
        double sin_lat         = std::sin(d2r * lat_y);
        double cos_lat         = std::cos(d2r * lat_y);
        double ZSYROT          = cos_south_pole_lat * sin_lat - sin_south_pole_lat * cos_lat * cos_lon_decr_sp;
        ZSYROT                 = std::max(std::min(ZSYROT, +1.0), -1.0);

        double PYROT = std::asin(ZSYROT) * r2d;

        double ZCYROT = std::cos(PYROT * d2r);
        double ZCXROT = (cos_south_pole_lat * cos_lat * cos_lon_decr_sp + sin_south_pole_lat * sin_lat) / ZCYROT;
        ZCXROT        = std::max(std::min(ZCXROT, +1.0), -1.0);
        double ZSXROT = cos_lat * sin_lon_decr_sp / ZCYROT;

        double PXROT = std::acos(ZCXROT) * r2d;

        if (ZSXROT < 0.0) {
            PXROT = -PXROT;
        }

        return PointLonLat(PXROT, PYROT);
    }

    PointLonLat magics_unrotate(const PointLonLat& point) const {
        double lat_y = point.lat();
        double lon_x = point.lon();

        double sin_south_pole_lat = std::sin(d2r * (south_pole_.lat() + 90.));
        double cos_south_pole_lat = std::cos(d2r * (south_pole_.lat() + 90.));
        double cos_lon            = std::cos(d2r * lon_x);
        double sin_lat            = std::sin(d2r * lat_y);
        double cos_lat            = std::cos(d2r * lat_y);
        double ZSYREG             = cos_south_pole_lat * sin_lat + sin_south_pole_lat * cos_lat * cos_lon;
        ZSYREG                    = std::max(std::min(ZSYREG, +1.0), -1.0);
        double PYREG              = std::asin(ZSYREG) * r2d;
        double ZCYREG             = std::cos(PYREG * d2r);
        double ZCXMXC             = (cos_south_pole_lat * cos_lat * cos_lon - sin_south_pole_lat * sin_lat) / ZCYREG;
        ZCXMXC                    = std::max(std::min(ZCXMXC, +1.0), -1.0);
        double ZSXMXC             = cos_lat * sin_lat / ZCYREG;
        double ZXMXC              = std::acos(ZCXMXC) * r2d;
        if (ZSXMXC < 0.0) {
            ZXMXC = -ZXMXC;
        }
        double PXREG = ZXMXC + south_pole_.lon();

        return PointLonLat(PXREG, PYREG);
    }
};

//-----------------------------------------------------------------------------


CASE("test_rotation_angle") {
    const int nx = 12, ny = 6;

    Rotation rot(Config("rotation_angle", 180.) | Config("north_pole", std::vector<double>{2.0, 46.7}));

    const PointLonLat ref[] = {
        {-178.00000,-46.70000}, {-178.00000,-76.70000},    {2.00000,-73.30000},    {2.00000,-43.30000},
           {2.00000,-13.30000},     {2.00000,16.70000},     {2.00000,46.70000}, {-178.00000,-46.70000},
         {140.11755,-68.00825},   {66.89106,-61.43199},   {40.42536,-36.43683},    {27.97634,-8.65459},
           {17.37657,19.46929},     {2.00000,46.70000}, {-178.00000,-46.70000},  {135.57504,-53.29508},
          {94.12083,-41.36507},   {69.20884,-20.05422},     {50.73654,3.83700},    {31.16557,27.31067},
            {2.00000,46.70000}, {-178.00000,-46.70000},  {141.90794,-39.07002},  {113.60146,-21.33906},
            {92.00000,0.00000},    {70.39854,21.33906},    {42.09206,39.07002},     {2.00000,46.70000},
        {-178.00000,-46.70000},  {152.83443,-27.31067},   {133.26346,-3.83700},   {114.79116,20.05422},
           {89.87917,41.36507},    {48.42496,53.29508},     {2.00000,46.70000}, {-178.00000,-46.70000},
         {166.62343,-19.46929},    {156.02366,8.65459},   {143.57464,36.43683},   {117.10894,61.43199},
           {43.88245,68.00825},     {2.00000,46.70000}, {-178.00000,-46.70000}, {-178.00000,-16.70000},
         {-178.00000,13.30000},  {-178.00000,43.30000},  {-178.00000,73.30000},     {2.00000,76.70000},
            {2.00000,46.70000}, {-178.00000,-46.70000}, {-162.62343,-19.46929},   {-152.02366,8.65459},
         {-139.57464,36.43683},  {-113.10894,61.43199},   {-39.88245,68.00825},     {2.00000,46.70000},
        {-178.00000,-46.70000}, {-148.83443,-27.31067},  {-129.26346,-3.83700},  {-110.79116,20.05422},
          {-85.87917,41.36507},   {-44.42496,53.29508},     {2.00000,46.70000}, {-178.00000,-46.70000},
        {-137.90794,-39.07002}, {-109.60146,-21.33906},    {-88.00000,0.00000},   {-66.39854,21.33906},
          {-38.09206,39.07002},     {2.00000,46.70000}, {-178.00000,-46.70000}, {-131.57504,-53.29508},
         {-90.12083,-41.36507},  {-65.20884,-20.05422},    {-46.73654,3.83700},   {-27.16557,27.31067},
            {2.00000,46.70000}, {-178.00000,-46.70000}, {-136.11755,-68.00825},  {-62.89106,-61.43199},
         {-36.42536,-36.43683},   {-23.97634,-8.65459},   {-13.37657,19.46929},     {2.00000,46.70000},
    };


    for (int i = 0, jglo = 0; i < nx; i++) {
        for (int j = 0; j < ny + 1; j++, jglo++) {
            double lon = static_cast<double>(i) * 360. / static_cast<double>(nx);
            double lat = static_cast<double>(j - ny / 2) * 90. / static_cast<double>(ny / 2);
            PointLonLat p0(lon, lat);
            PointLonLat p1 = rot.rotate(p0);
            PointLonLat p2 = rot.unrotate(p1);
            EXPECT(equivalent(p0, p2));
            EXPECT(equivalent(p1, ref[jglo]));
        }
    }
}

CASE("test_rotation_construction") {
    static const PointLonLat SP{0., -90.};
    static const PointLonLat NP{180., 90.};

    std::vector<PointLonLat> rotation_poles = {SP, NP, {0., -90.1}, {0., 90.1}};

    for (auto& p : rotation_poles) {
        Rotation s(Config("south_pole", std::vector<double>{p.lon(), p.lat()}));
        Log::info() << "rotate_south_pole=" << s << std::endl;
        EXPECT(s.rotated() == (p != SP));

        Rotation n(Config("north_pole", std::vector<double>{p.lon(), p.lat()}));
        Log::info() << "rotate_north_pole=" << n << std::endl;
        EXPECT(n.rotated() == (p != NP));
    }
}

CASE("test_rotation") {
    Config config;
    config.set("north_pole", std::vector<double>{-176, 40});
    Rotation rotation(config);
    MagicsRotation magics(rotation.southPole());
    Log::info() << rotation << std::endl;

    EXPECT(rotation.rotated());

    PointLonLat p, r;

    p = {0., 90.};
    r = {-176., 40.};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(magics.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
    EXPECT_EQUIVALENT(magics.unrotate(r), p);

    p = {0., 0.};
    r = {4., 50.};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(magics.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
    EXPECT_EQUIVALENT(magics.unrotate(r), p);

    p = {-180., 45.};
    r = {-176., -5.};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(magics.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
    EXPECT_EQUIVALENT(magics.unrotate(r), p);
}

CASE("test_no_rotation") {
    Config config;
    Rotation rotation(config);
    MagicsRotation magics(rotation.southPole());

    Log::info() << rotation << std::endl;

    EXPECT(not rotation.rotated());

    PointLonLat p, r;

    p = {0., 90.};
    r = p;
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(magics.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
    EXPECT_EQUIVALENT(magics.unrotate(r), p);

    p = {0., 0.};
    r = p;
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(magics.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
    EXPECT_EQUIVALENT(magics.unrotate(r), p);

    p = {-180., 45.};
    r = p;
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(magics.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
    EXPECT_EQUIVALENT(magics.unrotate(r), p);
}

CASE("test_rotation_angle_only") {
    Config config;
    config.set("rotation_angle", -180.);
    Rotation rotation(config);
    MagicsRotation magics(rotation.southPole());

    Log::info() << rotation << std::endl;

    EXPECT(rotation.rotated());

    PointLonLat p, r;

    p = {0., 90.};
    r = {-180., 90};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);

    p = {0., 0.};
    r = {-180., 0.};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);

    p = {270., 25.};
    r = {90., 25.};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);

    p = {-180., 45.};
    r = {-360., 45.};
    EXPECT_EQUIVALENT(rotation.rotate(p), r);
    EXPECT_EQUIVALENT(rotation.unrotate(r), p);
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
