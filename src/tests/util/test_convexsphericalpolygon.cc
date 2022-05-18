#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using ConvexSphericalPolygon = util::ConvexSphericalPolygon;

util::ConvexSphericalPolygon make_polygon(const std::initializer_list<PointLonLat>& list) {
    return util::ConvexSphericalPolygon{std::vector<PointLonLat>(list)};
}

template <class... Ts>
util::ConvexSphericalPolygon make_polygon(Ts&&... ts) {
    std::array<PointLonLat, sizeof...(Ts)> arr{ts...};
    return util::ConvexSphericalPolygon{arr};
}


util::ConvexSphericalPolygon make_polygon() {
    return util::ConvexSphericalPolygon{};
}

CASE("test default constructor") {
    ConvexSphericalPolygon p;
    EXPECT(bool(p) == false);
}

CASE("test_spherical_polygon_area") {
    auto plg1 = make_polygon({{0., 90.}, {0., 0.}, {90., 0.}});
    EXPECT_APPROX_EQ(plg1.area(), M_PI_2);
    auto plg2 = make_polygon({{0., 45.}, {0., 0.}, {90., 0.}, {90., 45.}});
    auto plg3 = make_polygon({{0., 90.}, {0., 45.}, {90., 45.}});
    Log::info() << "area diff: " << plg1.area() - plg2.area() - plg3.area() << std::endl;
    EXPECT_APPROX_EQ(std::abs(plg1.area() - plg2.area() - plg3.area()), 0, 1e-15);
}

CASE("test_spherical_polygon_intersection") {
    constexpr int nplg_f                             = 2;
    constexpr int nplg_g                             = 17;
    constexpr int nplg_i                             = nplg_f * nplg_g;
    std::array<ConvexSphericalPolygon, nplg_f> plg_f = {make_polygon({{0, 70}, {0, 60}, {40, 60}, {40, 70}}),
                                                        make_polygon({{0, 90}, {0, 0}, {40, 0}})};
    std::array<ConvexSphericalPolygon, nplg_g> plg_g = {
        make_polygon({{0, 60}, {0, 50}, {40, 60}}),  //0
        make_polygon({{0, 60}, {0, 50}, {20, 60}}),
        make_polygon({{10, 60}, {10, 50}, {30, 60}}),  //3
        make_polygon({{40, 80}, {0, 60}, {40, 60}}),
        make_polygon({{0, 80}, {0, 60}, {40, 60}}),  //5
        make_polygon({{20, 80}, {0, 60}, {40, 60}}),
        make_polygon({{20, 70}, {0, 50}, {40, 50}}),  //7
        make_polygon({{0, 90}, {0, 60}, {40, 60}}),
        make_polygon({{-10, 80}, {-10, 50}, {50, 80}}),  //9
        make_polygon({{0, 80}, {0, 50}, {40, 50}, {40, 80}}),
        make_polygon({{0, 65}, {20, 55}, {40, 60}, {20, 65}}),  //11
        make_polygon({{20, 65}, {0, 60}, {20, 55}, {40, 60}}),
        make_polygon({{10, 63}, {20, 55}, {30, 63}, {20, 65}}),  //13
        make_polygon({{20, 75}, {0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}}),
        make_polygon({{0, 50}, {0, 40}, {5, 45}}),  //15
        make_polygon({{0, 90}, {0, 80}, {20, 0}, {40, 80}}),
        make_polygon({{0, 65}, {0, 55}, {40, 65}, {40, 75}}),  //17
    };
    std::array<ConvexSphericalPolygon, nplg_i> plg_i = {
        make_polygon(),  //0
        make_polygon(),
        make_polygon(),  //2
        make_polygon({{0, 60}, {40, 60}, {40, 70}, {10, 70.8}}),
        make_polygon({{0, 70}, {0, 60}, {40, 60}, {30, 70.8}}),  //4
        make_polygon({{0, 60}, {40, 60}, {34.6, 70.5}, {5.3, 70.5}}),
        make_polygon({{7.5, 60.9}, {32.5, 60.9}, {20, 70}}),  //6
        make_polygon({{0, 70}, {0, 60}, {40, 60}, {40, 70}}),
        make_polygon({{0, 65.5}, {6.9, 70.6}, {0, 70}}),  //8
        make_polygon({{0, 70}, {0, 60}, {40, 60}, {40, 70}}),
        make_polygon({{0, 65}, {10, 61.1}, {40, 60}, {20, 65}}),  //10
        make_polygon({{0, 60}, {40, 60}, {20, 65}}),
        make_polygon({{12.6, 61.3}, {27.4, 61.3}, {30, 63}, {20, 65}, {10, 63}}),  //12
        make_polygon({{0, 70}, {1.9, 60.3}, {32.9, 60.9}, {40, 70}}),
        make_polygon(),  //14
        make_polygon({{13.7, 61.4}, {26.3, 61.4}, {30, 70.8}, {10, 70.8}}),
        make_polygon({{0, 65}, {0, 60}, {16.8, 61.5}, {40, 65}, {40, 70}, {15, 71.1}}),  //16
        make_polygon({{0, 60}, {0, 50}, {40, 60}}),
        make_polygon({{0, 60}, {0, 50}, {20, 60}}),  //18
        make_polygon({{10, 60}, {10, 50}, {30, 60}}),
        make_polygon({{0, 60}, {40, 60}, {40, 80}}),  //20
        make_polygon({{0, 80}, {0, 60}, {40, 60}}),
        make_polygon({{0, 60}, {40, 60}, {20, 80}}),  //22
        make_polygon({{0, 50}, {40, 50}, {20, 70}}),
        make_polygon({{0, 90}, {0, 60}, {40, 60}}),  //24
        make_polygon({{0, 65.5}, {40, 79.2}, {40, 80.8}, {0, 80.8}}),
        make_polygon({{0, 80}, {0, 50}, {40, 50}, {40, 80}}),  //26
        make_polygon({{0, 65}, {20, 55}, {40, 60}, {20, 65}}),
        make_polygon({{0, 60}, {20, 55}, {40, 60}, {20, 65}}),  //28
        make_polygon({{10, 63}, {20, 55}, {30, 63}, {20, 65}}),
        make_polygon({{0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}, {20, 75}}),  //30
        make_polygon({{0, 50}, {0, 40}, {5, 45}}),
        make_polygon({{0, 90}, {0, 80}, {20, 0}, {40, 80}}),  //32
        make_polygon({{0, 65}, {0, 55}, {40, 65}, {40, 75}})};
    for (int i = 0; i < nplg_f; i++) {
        for (int j = 0; j < nplg_g; j++) {
            Log::info() << "\n(" << i * nplg_g + j << ") Intersecting polygon\n    " << plg_f[i] << std::endl;
            Log::info() << "with polygon\n    " << plg_g[j] << std::endl;
            auto plg_fg = plg_f[i].intersect(plg_g[j]);
            auto plg_gf = plg_g[j].intersect(plg_f[i]);
            Log::info() << "got polygon\n    ";
            if (plg_fg) {
                Log::info() << plg_fg << std::endl;
                Log::info() << "    " << plg_gf << std::endl;
                EXPECT(plg_fg.equals(plg_gf));
                EXPECT(plg_fg.equals(plg_i[i * nplg_g + j], 0.1));
                Log::info() << "instead of polygon\n    ";
                Log::info() << plg_i[i * nplg_g + j] << std::endl;
            }
            else {
                Log::info() << "	empty" << std::endl;
                Log::info() << "instead of polygon\n    ";
                Log::info() << plg_i[i * nplg_g + j] << std::endl;
            }
        }
    }
}

CASE("Size of ConvexSphericalPolygon") {
    // This test illustrates that ConvexSphericalPolygon is allocated on the stack completely,
    // as sizeof(ConvexSphericalPolygon) includes space for MAX_SIZE coordinates of type PointXYZ
    EXPECT(sizeof(PointXYZ) == sizeof(double) * 3);
    size_t expected_size = 0;
    expected_size += (1 + ConvexSphericalPolygon::MAX_SIZE) * sizeof(PointXYZ);
    expected_size += sizeof(size_t);
    expected_size += sizeof(bool);
    expected_size += 2 * sizeof(double);
    EXPECT(sizeof(ConvexSphericalPolygon) >= expected_size);  // greater because compiler may add some padding
}

CASE("analyse intersect") {
    const double EPS  = std::numeric_limits<double>::epsilon();
    const double du   = 0.5;
    const double dv   = 1.1 * EPS;
    const double duc  = 0.5 * du;
    const double sduc = std::sqrt(1. - 0.25 * du * du);
    const double dvc  = 1. - 0.5 * dv * dv;
    const double sdvc = dv * std::sqrt(1. - 0.25 * dv * dv);
    PointXYZ s0p0{sduc, -duc, 0.};
    PointXYZ s0p1{sduc, duc, 0.};
    PointXYZ s1p0{dvc * sduc, -dvc * duc, -sdvc};
    PointXYZ s1p1{dvc * sduc, dvc * duc, sdvc};

    EXPECT_APPROX_EQ(dv, PointXYZ::norm(s0p0 - s1p0), EPS);
    EXPECT_APPROX_EQ(du, PointXYZ::norm(s0p0 - s0p1), EPS);
    EXPECT_APPROX_EQ(dv, PointXYZ::norm(s0p1 - s1p1), EPS);

    ConvexSphericalPolygon::GreatCircleSegment s1(s0p0, s0p1);
    ConvexSphericalPolygon::GreatCircleSegment s2(s1p0, s1p1);

    // analytical solution
    PointXYZ Isol{1., 0., 0.};

    // test "intersection"
    PointXYZ I = s1.intersect(s2);
    EXPECT_APPROX_EQ(std::abs(PointXYZ::norm(I) - 1.), 0., EPS);
    EXPECT_APPROX_EQ(PointXYZ::norm(I - Isol), 0., EPS);

    // test "contains"
    EXPECT(s1.contains(Isol) && s2.contains(Isol));
    EXPECT(s1.contains(I) && s2.contains(I));
}

CASE("source_covered") {
    const double dd = 0.;
    const auto csp0 = make_polygon({{0, 90}, {0, 0}, {90, 0}});
    double dcov1    = csp0.area();  // optimal coverage
    double dcov2    = dcov1;        // intersection-based coverage
    double dcov3    = dcov1;        // normalized intersection-based coverage
    double darea    = 0;            // commutative area error in intersection: |area(A^B)-area(B^A)|

    double max_tarea         = 0.;
    double min_tarea         = 0.;
    double accumulated_tarea = 0.;

    const int n       = 900;
    const int m       = 900;
    const double dlat = 90. / n;
    const double dlon = 90. / m;
    std::vector<ConvexSphericalPolygon> csp(n * m);
    std::vector<double> tgt_area(n * m);
    ATLAS_TRACE_SCOPE("1")
    for (int j = 0; j < m; j++) {
        csp[j] =
            make_polygon(PointLonLat{0, 90}, PointLonLat{dlon * j, 90 - dlat}, PointLonLat{dlon * (j + 1), 90 - dlat});
    }
    ATLAS_TRACE_SCOPE("2")
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < m; j++) {
            csp[i * m + j] = make_polygon(
                PointLonLat{dlon * j, 90 - dlat * i}, PointLonLat{dlon * j, 90 - dlat * (i + 1)},
                PointLonLat{dlon * (j + 1), 90 - dlat * (i + 1)}, PointLonLat{dlon * (j + 1), 90 - dlat * i});
        }
    }
    ConvexSphericalPolygon cspi0;
    ConvexSphericalPolygon csp0i;

    ATLAS_TRACE_SCOPE("3")
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            auto ipoly = i * m + j;
            cspi0      = csp[ipoly].intersect(csp0);  // intersect small csp[...] with large polygon csp0
                                                      // should be approx csp[...] as well
            csp0i = csp0.intersect(csp[ipoly]);       // opposite: intersect csp0 with small csp[...]
                                                      // should be approx csp[...] as well
            double a_csp   = csp[ipoly].area();
            double a_cspi0 = cspi0.area();  // should be approx csp[...].area()
            double a_csp0i = csp0i.area();  // should approx match a_cspi0
            EXPECT_APPROX_EQ(a_cspi0, a_csp, 1.e-10);
            EXPECT_APPROX_EQ(a_csp0i, a_csp, 1.e-10);
            darea = std::max(darea, a_cspi0 - a_csp0i);  // should remain approx zero
            dcov1 -= csp[ipoly].area();
            dcov2 -= a_cspi0;
            tgt_area[ipoly] = a_cspi0;
            accumulated_tarea += tgt_area[ipoly];
        }
    }

    // normalize weights
    double norm_fac = csp0.area() / accumulated_tarea;
    ATLAS_TRACE_SCOPE("4")
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            tgt_area[i * m + j] *= norm_fac;
            dcov3 -= tgt_area[i * m + j];
        }
    }

    Log::info() << " dlat, dlon : " << dlat << ", " << dlon << "\n";
    Log::info() << " max (commutative_area) : " << darea << "\n";
    Log::info() << " dcov1 : " << dcov1 << "\n";
    Log::info() << " dcov2 : " << dcov2 << "\n";
    Log::info() << " dcov3 : " << dcov3 << "\n";
    Log::info() << " accumulated small polygon area : " << accumulated_tarea << "\n";
    Log::info() << " large polygon area             : " << csp0.area() << "\n";
    EXPECT_APPROX_EQ(darea, 0., 1.e-10);
    EXPECT_APPROX_EQ(accumulated_tarea, csp0.area(), 1.e-8);
}

CASE("edge cases") {
    Log::info().precision(20);
    SECTION("CS-LFR-256 -> H1280 problem polygon intersection") {
        const auto plg0 = make_polygon({{-23.55468749999994, -41.11286269132660},
                                        {-23.20312500000000, -41.18816845938357},
                                        {-23.20312500000000, -40.83947225425061},
                                        {-23.55468749999994, -40.76429594967151}});
        const auto plg1 = make_polygon({{-23.30859375000000, -40.81704944888558},
                                        {-23.27343750000000, -40.85649237345376},
                                        {-23.23828125000000, -40.81704944888558},
                                        {-23.27343750000000, -40.77762996221442}});
        auto iplg       = plg0.intersect(plg1);
        auto jplg       = plg1.intersect(plg0);
        EXPECT_APPROX_EQ(iplg.area(), jplg.area(), 2e-12);  // can not take 1e-15
        EXPECT_EQ(iplg.size(), 3);
        EXPECT_EQ(jplg.size(), 3);
        EXPECT(iplg.equals(jplg, 5.e-7));
        Log::info() << "Intersection area difference: " << std::abs(iplg.area() - jplg.area()) << "\n";
    }
    SECTION("CS-LFR-16 -> O32 problem polygon intersection") {
        const auto plg0 = make_polygon({{174.3750000000001, -16.79832945594544},
                                        {174.3750000000001, -11.19720014633353},
                                        {168.7500000000000, -11.03919441545243},
                                        {168.7500000000000, -16.56868711281520}});
        const auto plg1 = make_polygon({{174.3750000000000, -12.55775611523100},
                                        {174.1935483870968, -15.34836475949100},
                                        {177.0967741935484, -15.34836475949100}});
        auto iplg       = plg0.intersect(plg1);
        auto jplg       = plg1.intersect(plg0);
        EXPECT_APPROX_EQ(iplg.area(), jplg.area(), 1e-13);  // can not take 1e-15
        EXPECT_EQ(iplg.size(), 3);
        EXPECT_EQ(jplg.size(), 3);
        EXPECT(iplg.equals(jplg, 1.e-11));
        Log::info() << "Intersection area difference: " << std::abs(iplg.area() - jplg.area()) << "\n";
    }
    SECTION("CS-LFR-2 -> O48 problem polygon intersection") {
        const auto plg0 = make_polygon({{180, -45}, {180, 0}, {135, 0}, {135, -35.26438968}});
        const auto plg1 = make_polygon({{180, -23.31573073}, {180, -25.18098558}, {-177.75, -23.31573073}});
        auto iplg       = plg0.intersect(plg1);
        auto jplg       = plg1.intersect(plg0);
        EXPECT_APPROX_EQ(iplg.area(), jplg.area(), 1e-15);
        EXPECT_EQ(iplg.size(), 0);
        EXPECT_EQ(jplg.size(), 0);
        EXPECT(iplg.equals(jplg, 1.e-12));
        Log::info() << "Intersection area difference: " << std::abs(iplg.area() - jplg.area()) << "\n";
    }
    SECTION("H128->H256 problem polygon intersection") {
        const auto plg0 = make_polygon({{42.45283019, 50.48004426},
                                        {43.33333333, 49.70239033},
                                        {44.15094340, 50.48004426},
                                        {43.26923077, 51.25580690}});
        const auto plg1 = make_polygon({{43.30188679, 50.48004426},
                                        {43.73831776, 50.09145680},
                                        {44.15094340, 50.48004426},
                                        {43.71428571, 50.86815893}});
        auto iplg       = plg0.intersect(plg1);
        auto jplg       = plg1.intersect(plg0);
        EXPECT_APPROX_EQ(iplg.area(), jplg.area(), 1e-15);
        EXPECT_EQ(iplg.size(), 4);
        EXPECT_EQ(jplg.size(), 4);
        EXPECT(iplg.equals(jplg, 1.e-12));
        Log::info() << "Intersection area difference: " << std::abs(iplg.area() - jplg.area()) << "\n";
    }

    SECTION("intesection of epsilon-distorted polygons") {
        const double eps = 1e-4;  // degrees
        const auto plg0  = make_polygon({{42.45283019, 50.48004426},
                                        {43.33333333, 49.70239033},
                                        {44.1509434, 50.48004426},
                                        {43.26923077, 51.2558069}});
        const auto plg1  = make_polygon({{42.45283019, 50.48004426 - eps},
                                        {43.33333333 + eps, 49.70239033},
                                        {44.1509434, 50.48004426 + eps},
                                        {43.26923077 - eps, 51.2558069}});
        auto iplg        = plg0.intersect(plg1);
        auto jplg        = plg1.intersect(plg0);
        EXPECT_APPROX_EQ(iplg.area(), jplg.area(), 1e-10);
        EXPECT_EQ(iplg.size(), 8);
        EXPECT_EQ(jplg.size(), 8);
        EXPECT(iplg.equals(jplg, 1.e-8));
        Log::info() << "Intersection area difference: " << std::abs(iplg.area() - jplg.area()) << "\n";
    }
}

//-----------------------------------------------------------------------------

}  // end namespace test
}  // end namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
