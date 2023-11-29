/*
 * (C) Copyright 2013 ECMWF
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */


#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

#include "atlas/grid.h"
#include "atlas/util/KDTree.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

constexpr double tol = 1.e-12;

//-----------------------------------------------------------------------------

CASE("Distributed creation from unstructured grid with halo") {
std::vector<PointXY> points = {
 {180,0},
 {90,0},
 {-90,0},
 {0,90},
 {0,-90},
 {0,0},
 {18,0},
 {36,0},
 {54,0},
 {72,0},
 {108,0},
 {126,0},
 {144,0},
 {162,0},
 {-162,0},
 {-144,0},
 {-126,0},
 {-108,0},
 {-72,0},
 {-54,0},
 {-36,0},
 {-18,0},
 {0,18},
 {0,36},
 {0,54},
 {0,72},
 {180,72},
 {180,54},
 {180,36},
 {180,18},
 {180,-18},
 {180,-36},
 {180,-54},
 {180,-72},
 {0,-72},
 {0,-54},
 {0,-36},
 {0,-18},
 {90,18},
 {90,36},
 {90,54},
 {90,72},
 {-90,72},
 {-90,54},
 {-90,36},
 {-90,18},
 {-90,-18},
 {-90,-36},
 {-90,-54},
 {-90,-72},
 {90,-72},
 {90,-54},
 {90,-36},
 {90,-18},
 {123.974,-58.6741},
 {154.087,-16.9547},
 {154.212,-58.8675},
 {114.377,-41.9617},
 {125.567,-23.5133},
 {137.627,-40.8524},
 {106.162,-24.5874},
 {158.508,-38.55},
 {137.826,-72.8109},
 {142.103,-26.799},
 {138.256,-13.8871},
 {168.39,-24.3266},
 {168.954,-12.0094},
 {117.333,-12.35},
 {102.254,-11.1537},
 {120.307,59.7167},
 {107.196,26.0167},
 {144.768,28.3721},
 {150.891,60.0343},
 {164.566,25.5053},
 {116.851,14.0295},
 {124.84,28.3978},
 {157.901,42.042},
 {111.41,43.1056},
 {134.333,44.6677},
 {103.277,11.707},
 {135.358,73.2119},
 {135.349,14.2311},
 {153.48,13.386},
 {168.071,11.5344},
 {-162.99,26.3775},
 {-147.519,56.1313},
 {-122.579,27.4824},
 {-117.909,59.2376},
 {-104.052,27.3616},
 {-153.107,14.9717},
 {-110.833,41.7436},
 {-144.847,32.8534},
 {-161.546,42.1031},
 {-129.866,44.5201},
 {-133.883,72.4163},
 {-166.729,11.8907},
 {-135.755,15.2529},
 {-106.063,14.4869},
 {-119.452,11.7037},
 {-146.026,-58.6741},
 {-115.913,-16.9547},
 {-115.788,-58.8675},
 {-155.623,-41.9617},
 {-144.433,-23.5133},
 {-132.373,-40.8524},
 {-163.838,-24.5874},
 {-111.492,-38.55},
 {-132.174,-72.8109},
 {-127.897,-26.799},
 {-131.744,-13.8871},
 {-101.61,-24.3266},
 {-101.046,-12.0094},
 {-152.667,-12.35},
 {-167.746,-11.1537},
 {-14.0127,-27.2963},
 {-59.193,-57.0815},
 {-56.465,-19.5751},
 {-27.056,-59.3077},
 {-57.124,-35.9752},
 {-33.4636,-28.3914},
 {-74.8037,-46.8602},
 {-40.089,-45.1376},
 {-74.8149,-28.3136},
 {-21.3072,-42.2177},
 {-44.0778,-72.6353},
 {-19.6969,-12.8527},
 {-40.1318,-12.1601},
 {-72.691,-11.4129},
 {-56.0261,58.6741},
 {-25.9127,16.9547},
 {-25.7876,58.8675},
 {-65.6229,41.9617},
 {-54.4335,23.5133},
 {-42.373,40.8524},
 {-73.838,24.5874},
 {-21.4917,38.55},
 {-42.1744,72.8109},
 {-37.8974,26.799},
 {-41.7437,13.8871},
 {-11.6095,24.3266},
 {-11.0459,12.0094},
 {-62.667,12.35},
 {-77.7456,11.1537},
 {30.3071,59.7167},
 {17.1956,26.0167},
 {54.7676,28.3721},
 {60.8915,60.0343},
 {74.5657,25.5053},
 {26.8506,14.0295},
 {34.8398,28.3978},
 {67.9014,42.042},
 {21.41,43.1056},
 {44.3335,44.6677},
 {13.2772,11.707},
 {45.3579,73.2119},
 {45.3492,14.2311},
 {63.4799,13.386},
 {78.0714,11.5344},
 {17.01,-26.3775},
 {32.4806,-56.1313},
 {57.4213,-27.4824},
 {62.0912,-59.2376},
 {75.9483,-27.3616},
 {26.893,-14.9717},
 {69.1672,-41.7436},
 {35.1527,-32.8534},
 {18.4543,-42.1031},
 {50.1344,-44.5201},
 {46.1172,-72.4163},
 {13.2711,-11.8907},
 {44.2448,-15.2529},
 {73.9368,-14.4869},
 {60.5478,-11.7037}
};
auto grid = UnstructuredGrid(points);

auto pointcloud = functionspace::PointCloud(grid,util::Config("halo_radius",2000000));

auto field = pointcloud.createField<double>(option::name("field"));

auto lonlat = array::make_view<double,2>(pointcloud.lonlat());
auto ghost  = array::make_view<int,1>(pointcloud.ghost());

auto view = array::make_view<double,1>(field);


auto fieldg_init = pointcloud.createField<double>(option::name("fieldg_init")|option::global());

if (mpi::rank() == 0) {
    auto viewg = array::make_view<double,1>(fieldg_init);
    gidx_t g=0;
    for (auto& p: grid.lonlat()) {
        double lat = p.lat() * M_PI/180.;
        viewg(g) = std::cos(4.*lat);
        g++;
    }
}

pointcloud.scatter(fieldg_init,field);

size_t count_ghost{0};
for (idx_t i=0; i<pointcloud.size(); ++i) {
    if( not ghost(i) ) {
        ++count_ghost;
        double lat = lonlat(i,1) * M_PI/180.;
        EXPECT_APPROX_EQ(view(i), std::cos(4.*lat), tol);
    }
    else {
        view(i) = 0.;
    }
}

if (mpi::size() == 4) {
    switch (mpi::rank()) {
        case 0: 
            EXPECT_EQ(pointcloud.size(),65);
            EXPECT_EQ(count_ghost,44);
            break;
        case 1:
            EXPECT_EQ(pointcloud.size(),69);
            EXPECT_EQ(count_ghost,43);
            break;
        case 2:
            EXPECT_EQ(pointcloud.size(),71);
            EXPECT_EQ(count_ghost,43);
            break;
        case 3:
            EXPECT_EQ(pointcloud.size(),61);
            EXPECT_EQ(count_ghost,43);
            break;
    }
}


field.haloExchange();

for (idx_t i=0; i<pointcloud.size(); ++i) {
    double lat = lonlat(i,1) * M_PI/180.;
    EXPECT_APPROX_EQ( view(i), std::cos(4.*lat), tol );
}



auto fieldg = pointcloud.createField<double>(option::name("field")|option::global());
if( mpi::rank() == 0 ) {
    EXPECT_EQ(fieldg.size(), grid.size());
}
else {
    EXPECT_EQ(fieldg.size(), 0);
}

pointcloud.gather(field, fieldg);

if (mpi::rank() == 0) {
    auto viewg = array::make_view<double,1>(fieldg);
    gidx_t g=0;
    for (auto& p: grid.lonlat()) {
        double lat = p.lat() * M_PI/180.;
        EXPECT_APPROX_EQ( viewg(g), std::cos(4.*lat), tol );
        g++;
    }
}

}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
