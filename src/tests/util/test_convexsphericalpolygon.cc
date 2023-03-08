#include <array>
#include <cmath>
#include <iostream>
#include <vector>
#include <type_traits>
#include <initializer_list>

#include "atlas/util/ConvexSphericalPolygon.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/Point.h"

#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {

using ConvexSphericalPolygon = util::ConvexSphericalPolygon;
const double EPS = std::numeric_limits<double>::epsilon();

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


template<typename It>
std::string to_json(const It& begin, const It& end, int precision = 0) {
    std::stringstream ss;
    ss << "[\n";
    size_t size = std::distance(begin,end);
    size_t c=0;
    for( auto it = begin; it != end; ++it, ++c ) {
        ss << "  " << it->json(precision);
        if( c < size-1 ) {
            ss << ",\n";
        }
    }
    ss << "\n]";
    return ss.str();
}

template<typename ConvexSphericalPolygonContainer>
std::string to_json(const ConvexSphericalPolygonContainer& polygons, int precision = 0) {
    return to_json(polygons.begin(),polygons.end(),precision);
}
std::string to_json(std::initializer_list<ConvexSphericalPolygon>&& polygons, int precision = 0) {
    return to_json(polygons.begin(),polygons.end(),precision);
}

void check_intersection(const ConvexSphericalPolygon& plg1, const ConvexSphericalPolygon& plg2, const ConvexSphericalPolygon& iplg_sol, double pointsSameEPS = 5.e6 * EPS, std::ostream* out = nullptr) {
    auto iplg       = plg1.intersect(plg2, out, pointsSameEPS);
    Log::info().indent();
    Log::info() << "plg1 area : " << plg1.area() << "\n";
    Log::info() << "plg2 area : " << plg2.area() << "\n";
    Log::info() << "iplg area : " << iplg.area() << "\n";
    Log::info() << "json: \n" << to_json({plg1,plg2,iplg},20) << "\n";
    EXPECT(std::min(plg1.area(), plg2.area()) >= iplg.area());
    EXPECT(iplg.equals(iplg_sol, 1.e-8));
    Log::info().unindent();
}

CASE("test default constructor") {
    ConvexSphericalPolygon p;
    EXPECT(bool(p) == false);
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
    PointXYZ I12 = s1.intersect(s2);
    PointXYZ I21 = s2.intersect(s1);
    EXPECT_APPROX_EQ(std::abs(PointXYZ::norm(I12) - 1.), 0., EPS);
    EXPECT_APPROX_EQ(PointXYZ::norm(I12 - Isol), 0., EPS);
    EXPECT_APPROX_EQ(std::abs(PointXYZ::norm(I21) - 1.), 0., EPS);
    EXPECT_APPROX_EQ(PointXYZ::norm(I21 - Isol), 0., EPS);
}

CASE("test_json_format") {
    auto plg = make_polygon({{0., 45.}, {0., 0.}, {90., 0.}, {90., 45.}});
    EXPECT_EQ(plg.json(), "[[0,45],[0,0],[90,0],[90,45]]");

    std::vector<ConvexSphericalPolygon> plg_g = {
        make_polygon({{0, 60}, {0, 50}, {40, 60}}),  //0
        make_polygon({{0, 60}, {0, 50}, {20, 60}}),
        make_polygon({{10, 60}, {10, 50}, {30, 60}}),  //3
    };
    Log::info() << to_json(plg_g) << std::endl;

}

CASE("test_spherical_polygon_area") {
    auto plg1 = make_polygon({{0., 90.}, {0., 0.}, {90., 0.}});
    EXPECT_APPROX_EQ(plg1.area(), M_PI_2);
    auto plg2 = make_polygon({{0., 45.}, {0., 0.}, {90., 0.}, {90., 45.}});
    auto plg3 = make_polygon({{0., 90.}, {0., 45.}, {90., 45.}});
    Log::info().indent();
    EXPECT(plg1.area() - plg2.area() - plg3.area() <= EPS);
    Log::info().unindent();
    EXPECT_APPROX_EQ(std::abs(plg1.area() - plg2.area() - plg3.area()), 0, 1e-15);
}

CASE("test_spherical_polygon_intersection") {
    constexpr int nplg_f                             = 2;
    constexpr int nplg_g                             = 17;
    constexpr int nplg_i                             = nplg_f * nplg_g;
    std::array<ConvexSphericalPolygon, nplg_f> plg_f = {make_polygon({{0, 70}, {0, 60}, {40, 60}, {40, 70}}),
                                                        make_polygon({{0, 90}, {0, 0}, {40, 0}})};
    std::array<ConvexSphericalPolygon, nplg_g> plg_g = {
        make_polygon({{0, 60}, {0, 50}, {40, 60}}),                            // 0
        make_polygon({{0, 60}, {0, 50}, {20, 60}}),                            // 1
        make_polygon({{10, 60}, {10, 50}, {30, 60}}),                          // 2
        make_polygon({{40, 80}, {0, 60}, {40, 60}}),                           // 3
        make_polygon({{0, 80}, {0, 60}, {40, 60}}),                            // 4
        make_polygon({{20, 80}, {0, 60}, {40, 60}}),                           // 5
        make_polygon({{20, 70}, {0, 50}, {40, 50}}),                           // 6
        make_polygon({{0, 90}, {0, 60}, {40, 60}}),                            // 7
        make_polygon({{-10, 80}, {-10, 50}, {50, 80}}),                        // 8
        make_polygon({{0, 80}, {0, 50}, {40, 50}, {40, 80}}),                  // 9
        make_polygon({{0, 65}, {20, 55}, {40, 60}, {20, 65}}),                 // 10
        make_polygon({{20, 65}, {0, 60}, {20, 55}, {40, 60}}),                 // 11
        make_polygon({{10, 63}, {20, 55}, {30, 63}, {20, 65}}),                // 12
        make_polygon({{20, 75}, {0, 70}, {5, 5}, {10, 0}, {20, 0}, {40, 70}}), // 13
        make_polygon({{0, 50}, {0, 40}, {5, 45}}),                             // 14
        make_polygon({{0, 90}, {0, 80}, {20, 0}, {40, 80}}),                   // 15
        make_polygon({{0, 65}, {0, 55}, {40, 65}, {40, 75}}),                  // 16
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
    Log::info().indent();
    for (int i = 0; i < nplg_f; i++) {
        for (int j = 0; j < nplg_g; j++) {
            auto plg_fg = plg_f[i].intersect(plg_g[j], nullptr, 5.e6 * std::numeric_limits<double>::epsilon());
            bool polygons_equal = plg_i[i * nplg_g + j].equals(plg_fg, 0.1);
            EXPECT(polygons_equal);
            if( not polygons_equal or (i==0 && j==10)) {
                Log::info() << "\nIntersected the polygon plg_f["<<i<<"]\n\t" << plg_f[i].json() << std::endl;
                Log::info() << "with the polygon plg_g["<<j<<"]\n\t" << plg_g[j].json() << std::endl;
                Log::info() << "and got the polygon plg_fg\n\t";
                Log::info() << plg_fg.json() << std::endl;
                Log::info() << "instead of the polygon plg_i["<<i*nplg_g+j<<"] \n\t";
                Log::info() << plg_i[i * nplg_g + j].json() << std::endl;
                Log::info() << "polygons.json ( [ plg_f[i], plg_g[j], plg_fg, plg_i["<<i*nplg_g+j<<"] ] ) = \n" << to_json({plg_f[i],plg_g[j],plg_fg,plg_i[i*nplg_g+j]}) << std::endl;
            }
        }
    }
    Log::info().unindent();
}

CASE("test_spherical_polygon_intersection_stretched") {
    // plg1 is an octant
    const auto plg1 = make_polygon({{0, 90}, {0, 0}, {90, 0}});
    Log::info().precision(20);

    const double delta = 1.1 * EPS;
    const double u = 1. - delta;
    const double v = std::sqrt(1 - u * u);
    const double w = delta;
    const double a = sqrt(1 - w * w);

    // plg2 is a stretched polygon
    std::vector<PointXYZ> plg2_points = {
        PointXYZ{u,   v,  0.},
        PointXYZ{a,   0., w },
        PointXYZ{u,  -v,  0.},
        PointXYZ{0., 0., -1.}};
    auto plg2 = util::ConvexSphericalPolygon(plg2_points.data(), plg2_points.size());
    EXPECT(plg2.size() == 4);

    auto plg11 = plg1.intersect(plg1);
    auto plg12 = plg1.intersect(plg2);
    auto plg22 = plg2.intersect(plg2);

    EXPECT(plg11.size() == 3);
    EXPECT(plg12.size() == 3);
    EXPECT(plg22.size() == 4);

    // plg3 is the analytical solution
    std::vector<PointXYZ> plg3_points = {
        PointXYZ{u,  v, 0.},
        PointXYZ{a,  0, w},
        PointXYZ{1, 0, 0.}};
    auto plg3 = util::ConvexSphericalPolygon(plg3_points.data(), plg3_points.size());
    EXPECT(plg3.size() == 3);
    EXPECT(plg12.size() == 3);
    auto plg_sol = make_polygon({{1.2074e-06, 0.}, {0., 1.3994e-14}, {0., 0.}});
    EXPECT(plg12.equals(plg_sol, 1e-10));

    Log::info().indent();
    Log::info() << "delta                       : " << delta << std::endl;
    Log::info() << "plg12.area                  : " << plg12.area() << std::endl;
    Log::info() << "exact intersection area     : " << plg3.area() << std::endl;
    double error_area = plg12.area() - plg3.area();
    EXPECT(error_area < EPS and error_area >= 0.);
    Log::info().unindent();
    Log::info().precision(-1);
}

CASE("test_lonlat_pole_problem") {
    const auto north_octant = make_polygon({{0, 90}, {0, 0}, {90, 0}});
    const double first_lat = 90. - 1.e+12 * EPS;
    const int m = 10000;
    const double dlon = 90. / m;
    std::vector<ConvexSphericalPolygon> csp(m);
    Log::info().indent();

    ATLAS_TRACE_SCOPE("create polygons")
    for (int j = 0; j < m; j++) {
        csp[j] =
            make_polygon(PointLonLat{0, 90}, PointLonLat{dlon * j, first_lat}, PointLonLat{dlon * (j + 1), first_lat});
    }

    double max_area_overshoot = 0.;
    int false_zero = 0;
    for (int j = 0; j < m; j++) {
        auto csp_i = csp[j].intersect(north_octant);
        // intersection area should not be larger than its father polygon's
        max_area_overshoot = std::max(max_area_overshoot, csp_i.area() - csp[j].area());
        if (csp_i.area() < EPS) {
            false_zero++;
        }
    }
    Log::info() << "False zero area intersection: " << false_zero << std::endl;
    Log::info() << "Max area overshoot: " << max_area_overshoot << std::endl;
    EXPECT(max_area_overshoot <= m * EPS);
    EXPECT(false_zero == 0);
    Log::info().unindent();
}

CASE("test_thin_elements_area") {
    const auto north_octant = make_polygon({{0, 90}, {0, 0}, {90, 0}});
    const auto south_octant = make_polygon({{0,0}, {0, -90},{90, 0}});
    const int n       = 3;
    const int m       = 2500;
    const double dlat = 180. / n;
    const double dlon = 90. / m;

    ATLAS_ASSERT(n > 1);
    std::vector<ConvexSphericalPolygon> csp(n * m);

    ATLAS_TRACE_SCOPE("create polygons")
    for (int j = 0; j < m; j++) {
        csp[j] =
            make_polygon(PointLonLat{0, 90}, PointLonLat{dlon * j, 90 - dlat}, PointLonLat{dlon * (j + 1), 90 - dlat});
    }
    for (int i = 1; i < n - 1; i++) {
        for (int j = 0; j < m; j++) {
            csp[i * m + j] = make_polygon(
                PointLonLat{dlon * j, 90 - dlat * i}, PointLonLat{dlon * j, 90 - dlat * (i + 1)},
                PointLonLat{dlon * (j + 1), 90 - dlat * (i + 1)}, PointLonLat{dlon * (j + 1), 90 - dlat * i});
        }
    }
    for (int j = 0; j < m; j++) {
        csp[(n - 1) * m + j] =
            make_polygon(PointLonLat{dlon * j, 90 - dlat * (n-1)}, PointLonLat{dlon * j, -90.},
                PointLonLat{dlon * (j + 1), 90 - dlat * (n-1)});
    }

    double coverage_north       = 0.;   // north octant coverage by intersections with "csp"s
    double coverage_south       = 0.;   // south octant coverage by intersections with "csp"s
    double coverage_norm        = 0.;   // area sum of intersections in the north octant normalised to sum up to the area of the north octant
    double coverage_csp         = 0.;   // area sum of all "csp"s
    double accumulated_tarea    = 0.;
    std::vector<double> i_north_area(n * m);

    ATLAS_TRACE_SCOPE("intersect polygons")
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            auto ipoly = i * m + j;
            double a_north_csp_i    = csp[ipoly].intersect(north_octant).area();  // intersect narrow csp with the north octant
            double a_south_csp_i    = csp[ipoly].intersect(south_octant).area();  // intersect narrow csp with the south octant
            if (i == 0) {
                if (n == 2) {
                    EXPECT_APPROX_EQ(a_north_csp_i, csp[ipoly].area(), EPS);
                }
                else {
                    if (a_north_csp_i > csp[ipoly].area()) {
                        Log::info() << " error: " << a_north_csp_i - csp[ipoly].area() << std::endl;
                    }
                }
            }
            coverage_north      +=  a_north_csp_i;
            coverage_csp        +=  csp[ipoly].area();
            coverage_south      +=  a_south_csp_i;
            i_north_area[ipoly] =   a_north_csp_i;
        }
    }

    // normalise weights of the intersection polygons to sum up to the area of the north octant
    double norm_fac = north_octant.area() / coverage_north;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            coverage_norm += i_north_area[i * m + j] * norm_fac;
        }
    }
    EXPECT(north_octant.area() - M_PI_2 < EPS); // spherical area error for the north octant
    EXPECT(south_octant.area() - M_PI_2 < EPS); // spherical area error for the south octant
    EXPECT(north_octant.area() - coverage_north < m * EPS); // error polygon intersec. north
    EXPECT(south_octant.area() - coverage_south < m * EPS); // error polygon intersec. north
    EXPECT(north_octant.area() - coverage_norm < m * EPS); // error polygon intersec. north
    auto err = std::max(north_octant.area() - coverage_north, south_octant.area() - coverage_south);
    auto err_normed = north_octant.area() - coverage_norm;
    Log::info() << "Octant coverting error        : " << err << std::endl;
    Log::info() << "Octant coverting error normed : " << err_normed << std::endl;
}

CASE("intesection of epsilon-distorted polygons") {
        const double eps = 1e-4;  // degrees
        const auto plg0     = make_polygon({{42.45283019, 50.48004426},
                                            {43.33333333, 49.70239033},
                                            {44.1509434, 50.48004426},
                                            {43.26923077, 51.2558069}});
        const auto plg1     = make_polygon({{42.45283019, 50.48004426 - eps},
                                            {43.33333333 + eps, 49.70239033},
                                            {44.1509434, 50.48004426 + eps},
                                            {43.26923077 - eps, 51.2558069}});
        const auto iplg_sol = make_polygon({{44.15088878324276, 50.48009332686897},
                                            {43.68455392823953, 50.89443301919586},
                                            {43.26918271448949, 51.25576215711414},
                                            {42.86876285000331, 50.87921047438197},
                                            {42.45288468219661, 50.47999711267543},
                                            {42.92307395320301, 50.06869923211562},
                                            {43.33338148668725, 49.70243705225555},
                                            {43.72937844034824, 50.08295071539503}});
        check_intersection(plg0, plg1, iplg_sol);
}

CASE("Edge cases in polygon intersections") {
    Log::info().precision(20);

    SECTION("CS-LFR-256 -> H1280") {
        const auto plg0 = make_polygon({{-23.55468749999994, -41.11286269132660},
                                        {-23.20312500000000, -41.18816845938357},
                                        {-23.20312500000000, -40.83947225425061},
                                        {-23.55468749999994, -40.76429594967151}});
        const auto plg1 = make_polygon({{-23.30859375000000, -40.81704944888558},
                                        {-23.27343750000000, -40.85649237345376},
                                        {-23.23828125000000, -40.81704944888558},
                                        {-23.27343750000000, -40.77762996221442}});
        auto iplg       = plg0.intersect(plg1);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg},16) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        Log::info().unindent();
    }
    SECTION("CS-LFR-16 -> O32") {
        const auto plg0 = make_polygon({{174.3750000000001, -16.79832945594544},
                                        {174.3750000000001, -11.19720014633353},
                                        {168.7500000000000, -11.03919441545243},
                                        {168.7500000000000, -16.56868711281520}});
        const auto plg1 = make_polygon({{174.3750000000000, -12.55775611523100},
                                        {174.1935483870968, -15.34836475949100},
                                        {177.0967741935484, -15.34836475949100}});
        auto iplg       = plg0.intersect(plg1);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg},16) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        Log::info().unindent();
    }
    SECTION("CS-LFR-2 -> O48") {
        const auto plg0 = make_polygon({{180, -45}, {180, 0}, {135, 0}, {135, -35.26438968}});
        const auto plg1 = make_polygon({{180, -23.31573073}, {180, -25.18098558}, {-177.75, -23.31573073}});
        auto iplg       = plg0.intersect(plg1);
        EXPECT_EQ(iplg.size(), 0);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg},10) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        EXPECT_APPROX_EQ(iplg.area(), 0.);
        Log::info().unindent();
    }
    SECTION("H128->H256") {
        const auto plg0 = make_polygon({{42.45283019, 50.48004426},
                                        {43.33333333, 49.70239033},
                                        {44.15094340, 50.48004426},
                                        {43.26923077, 51.25580690}});
        const auto plg1 = make_polygon({{43.30188679, 50.48004426},
                                        {43.73831776, 50.09145680},
                                        {44.15094340, 50.48004426},
                                        {43.71428571, 50.86815893}});
        auto iplg       = plg0.intersect(plg1);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg},10) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        Log::info().unindent();
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
        EXPECT_EQ(iplg.size(), 8);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg},10) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        Log::info().unindent();
    }
    SECTION("one") { // TODO: remove, do not know what example
        auto plg0 = make_polygon({{0,90},{67.463999999999998636,89.999899999085130275},{67.5,89.999899999085130275}});
        auto plg1 = make_polygon({{72,85.760587120443801723},{90,85.760587120443801723},{-90,85.760587120443801723}});
        auto iplg       = plg0.intersect(plg1);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg},20) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        Log::info().unindent();
    }

    SECTION("debug") { // TODO: remove, do not know what example
        auto plg0 = make_polygon({{168.599999999999994,-6.88524379355500038},{168.561872909699019,-7.16627415158899961},{168.862876254180634,-7.16627415158899961}});
        auto plg1 = make_polygon({{168.84,-6.84},{168.84,-7.2},{169.2,-7.2},{169.2,-6.84}});
        auto plg2 = make_polygon({{168.48,-6.84},{168.48,-7.2},{168.84,-7.2},{168.84,-6.84}});
        auto iplg       = plg0.intersect(plg1);
        auto iplg2       = plg0.intersect(plg2);
        Log::info().indent();
        Log::info() << "source area: " << plg0.area() << "\n";
        Log::info() << "target area: " << plg1.area() << "\n";
        Log::info() << "inters area: " << iplg.area() << "\n";
        Log::info() << "json: \n" << to_json({plg0,plg1,iplg2,iplg,iplg2},20) << "\n";
        EXPECT(std::min(plg0.area(), plg1.area()) >= iplg.area());
        Log::info().unindent();
    }


    SECTION("debug2") { // TODO: remove, do not know what example
        auto plg0           = make_polygon({{168.599999999999994, -6.88524379355500038},
                                            {168.561872909699019, -7.16627415158899961},
                                            {168.862876254180634, -7.16627415158899961}});
        auto plg1           = make_polygon({{168.48, -6.84},
                                            {168.48, -7.20},
                                            {168.84, -7.20},
                                            {168.84, -6.84}});
        const auto iplg_sol = make_polygon({{168.600000000000, -6.885243793555000},
                                            {168.561872909699, -7.166274151589000},
                                            {168.840000000000, -7.166281023991750},
                                            {168.840000000000, -7.141837487873564}});
        check_intersection(plg0, plg1, iplg_sol);
    }

    SECTION("O320c_O320n_tcell-2524029") {
        auto plg0           = make_polygon({{ 16.497973615369710, 89.85157892074884},
                                            {  0.000000000000000, 89.87355342974176},
                                            {-54.000000000000010, 89.78487690721863},
                                            { -9.000000000000002, 89.84788464490568}});
        auto plg1           = make_polygon({{ 36.000000000000000, 89.78487690721863},
                                            {-54.000000000000010, 89.78487690721863},
                                            {-36.000000000000000, 89.78487690721863}});
        const auto iplg_sol = make_polygon();
        check_intersection(plg0, plg1, iplg_sol);
    }

    SECTION("O128c_F128c_tcell-77536") {
        auto plg0           = make_polygon({{157.5,-16.49119584364},{157.5,-17.192948774143},{158.203125,-17.192948774143},{158.203125,-16.49119584364}});
        auto plg1           = make_polygon({{157.7064220183486,-16.49119584364},{157.5,-17.192948774143},{158.3333333333333,-17.192948774143}});
        const auto iplg_sol = make_polygon({{157.7066386314724, -16.49143953797217},
                                            {157.7063506671565, -16.49143933951490},
                                            {157.5000000000000, -17.19294877414300},
                                            {158.2031250000000, -17.19294877414292},
                                            {158.2031250000000, -17.04778221576889}});
        check_intersection(plg0, plg1, iplg_sol);
    }

    SECTION("O128c_F128c_tcell") {
        auto plg0           = make_polygon({{135,-10.505756145244},{135,-11.906523334954},{136.40625,-11.906523334954},{136.40625,-10.505756145244}});
        auto plg1           = make_polygon({{135.7377049180328,-10.505756145244},{135,-11.906523334954},{136.5,-11.906523334954}});
        const auto iplg_sol = make_polygon({{135.7381225488238, -10.50652773728132},
                                            {135.7373006953013, -10.50652782622945},
                                            {135.0000000000000, -11.90652333495400},
                                            {136.4062500000000, -11.90652333495396},
                                            {136.4062500000000, -11.73509894855489}});
        check_intersection(plg0, plg1, iplg_sol);
    }

    SECTION("O128c_F128c_tcell-2") {
        auto plg0           = make_polygon({{134.296875,-10.877172064989},{134.296875,-11.578925065131},{135,-11.578925065131},{135,-10.877172064989}});
        auto plg1           = make_polygon({{134.6153846153846,-10.877172064989},{134.2241379310345,-11.578925065131},{135,-11.578925065131}});
        const auto iplg_sol = make_polygon({{134.6154929040914, -10.87737018871372},
                                            {134.6152744739456, -10.87737016536156},
                                            {134.2968750000000, -11.44875997313061},
                                            {134.2968749999999, -11.57892506513095},
                                            {135.0000000000000, -11.57892506513100}});
        check_intersection(plg0, plg1, iplg_sol);
    }

}

//-----------------------------------------------------------------------------

}  // end namespace test
}  // end namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
