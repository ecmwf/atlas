#include "atlas/grid/CubedSphereGrid2.h"
#include "atlas/grid/Iterator.h"
#include "eckit/geometry/Point2.h"
#include "tests/AtlasTestEnvironment.h"

namespace atlas {
namespace test {
namespace {
using Point2 = eckit::geometry::Point2;

template <typename point2_derived>
bool compare_2D_points(std::vector<point2_derived> a, std::vector<Point2> b, double tolerance = 1e-4) { 
    // Uses a tolerance as the values are stored more precisely than they are printed
    ATLAS_ASSERT(a.size() == b.size());
    bool equal = true;
    for (int i = 0; i < b.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            if (std::abs(b[i][j] - a[i][j]) > tolerance) {
                std::cout << "[" << i << ", " << j << "]\n\t" << a[i][j] << " != " << b[i][j]
                          << "\n\tdiff = " << b[i][j] - a[i][j] << std::endl;
                equal = false;
            }
        }
    }
    return equal;
}

CASE("cubed_sphere_instantiation") {
    const int n = 2;
    const Grid grid = CubedSphereGrid2(n);

    EXPECT(grid.name() == "CS-LFR-" + std::to_string(n) + "-2");
    EXPECT(grid.type() == "cubedsphere2");
    EXPECT(grid.size() == n * n * 6);
}

CASE("cubed_sphere_grid_kgo") {
    // Lonlat and XY are both currently lonlat positions
    std::vector<Point2> kgo_lonlat { // N = 2
        {-22.5,20.941},  {22.5,20.941},   {-22.5,-20.941},  {22.5,-20.941},  {67.5,20.941},   {112.5,20.941},
        {67.5,-20.941},  {112.5,-20.941}, {157.5,20.941},   {-157.5,20.941}, {157.5,-20.941}, {-157.5,-20.941},
        {-112.5,20.941}, {-67.5,20.941},  {-112.5,-20.941}, {-67.5,-20.941}, {-45,59.6388},   {-135,59.6388},
        {45,59.6388},    {135,59.6388},   {45,-59.6388},    {135,-59.6388},  {-45,-59.6388},  {-135,-59.6388}
    };

    const Grid grid = CubedSphereGrid2(2);

    // LonLat
    std::cout << "Checking lonlat" << std::endl;
    std::vector<PointLonLat> points_lonlat;
    for (const auto &ll : grid.lonlat()) {
        points_lonlat.push_back(ll);
    }
    EXPECT(points_lonlat.size() == static_cast<size_t>(grid.size()));
    EXPECT(compare_2D_points<PointLonLat>(points_lonlat, kgo_lonlat));

    // XY
    std::cout << "Checking xy" << std::endl;
    std::vector<PointXY> points_xy;
    for (const auto &xy : grid.xy()) {
        points_xy.push_back(xy);
    }
    EXPECT(points_xy.size() == static_cast<size_t>(grid.size()));
    EXPECT(compare_2D_points<PointXY>(points_xy, kgo_lonlat));
}

}  // namespace
}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
