/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "atlas/grid.h"
#include "atlas/util/KDTree.h"

#include "tests/AtlasTestEnvironment.h"

using namespace atlas::util;

namespace atlas {
namespace test {

//------------------------------------------------------------------------------------------------
namespace {  // helpers

template <typename T>
std::string to_string(const std::vector<T>& vector) {
    std::stringstream s;
    s << vector;
    return s.str();
}
template <typename T>
std::string to_string(const T& v) {
    return std::to_string(v);
}

class PayloadGenerator {
public:
    class iterator {
        friend class PayloadGenerator;

    public:
        long int operator*() const { return i_; }
        const iterator& operator++() {
            ++i_;
            return *this;
        }
        /*
        // unused
        iterator operator++( int ) {
            iterator copy( *this );
            ++i_;
            return copy;
        }*/

        bool operator==(const iterator& other) const { return i_ == other.i_; }
        bool operator!=(const iterator& other) const { return i_ != other.i_; }

    protected:
        iterator(long int start): i_(start) {}

    private:
        unsigned long i_;
    };

    iterator begin() const { return begin_; }
    iterator end() const { return end_; }
    PayloadGenerator(long int end): begin_(0), end_(end) {}

    template <typename Container, typename Value = typename Container::value_type>
    void fill(Container& container) {
        std::iota(container.begin(), container.end(), Value(*begin_));
    }

    template <typename Value, size_t size>
    std::array<Value, size> make_array() {
        std::array<Value, size> array;
        fill(array);
        return array;
    }

private:
    iterator begin_;
    iterator end_;
};

static double radius() {
    return util::Earth::radius();
}

static Geometry& geometry() {
    static Geometry _geometry(radius());
    return _geometry;
}

static PointXYZ make_xyz(const PointLonLat& lonlat) {
    return geometry().xyz(lonlat);
}

static std::array<double, 7>& test_lon() {
    static std::array<double, 7> lon = {0., 30., 60., 90., 120., 150., 180.};
    return lon;
}

static std::array<double, 7>& test_lat() {
    static std::array<double, 7> lat = {90., 60., 30., 0., -30., -60., -90.};
    return lat;
}

static std::array<PointLonLat, 7>& test_lonlat() {
    static std::array<PointLonLat, 7> lonlat{PointLonLat{0., 90.},   PointLonLat{30., 60.},   PointLonLat{60., 30.},
                                             PointLonLat{90., 0.},   PointLonLat{120., -30.}, PointLonLat{150., -60.},
                                             PointLonLat{180., -90.}};
    return lonlat;
}

static std::array<PointXYZ, 7>& test_xyz() {
    static std::array<PointXYZ, 7> xyz{make_xyz(PointLonLat{0., 90.}),    make_xyz(PointLonLat{30., 60.}),
                                       make_xyz(PointLonLat{60., 30.}),   make_xyz(PointLonLat{90., 0.}),
                                       make_xyz(PointLonLat{120., -30.}), make_xyz(PointLonLat{150., -60.}),
                                       make_xyz(PointLonLat{180., -90.})};
    return xyz;
}

static std::array<idx_t, 7>& test_payloads() {
    static auto payloads = PayloadGenerator(7).make_array<idx_t, 7>();
    EXPECT_EQ(std::distance(payloads.begin(), payloads.end()), 7);
    return payloads;
}

void validate(const KDTree<idx_t>& tree) {
    EXPECT_NO_THROW(tree.closestPoint(PointLonLat{180., 45.}));
    // Search 4 nearest neighbours (k=4), sorted by shortest distance
    auto neighbours        = tree.closestPoints(PointLonLat{89.9, 44.9}, 4);
    auto expected_payloads = std::vector<idx_t>{2, 1, 3, 0};
    EXPECT_EQ(neighbours.payloads(), expected_payloads);
}

static const IndexKDTree& search() {
    static IndexKDTree kdtree = []() {
        IndexKDTree kdtree(geometry());
        auto grid = Grid{"O32"};
        kdtree.build(grid.lonlat(), PayloadGenerator(grid.size()));
        return kdtree;
    }();
    return kdtree;
}

}  // namespace
//------------------------------------------------------------------------------------------------

static bool ECKIT_515_implemented = ATLAS_ECKIT_VERSION_AT_LEAST(1, 13, 2);
// --> implements eckit::KDTree::size() and eckit::KDTree::empty()

CASE("test kdtree") {
    auto grid = Grid{"O32"};

    IndexKDTree search(geometry());
    EXPECT(search.empty());

    search.reserve(grid.size());
    idx_t n{0};
    for (auto& point : grid.lonlat()) {
        search.insert(point, n++);
        if (ECKIT_515_implemented) {
            EXPECT(search.empty());
        }
    }
    search.build();
    EXPECT_EQ(search.size(), grid.size());
    EXPECT_NO_THROW(search.closestPoint(PointLonLat{180., 45.}));
    // ...
    // Search 4 nearest neighbours (k=4), sorted by shortest distance
    auto neighbours          = search.closestPoints(PointLonLat{180., 45.}, 4).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761};
    EXPECT_EQ(neighbours, expected_neighbours);
}


CASE("test assertion") {
    auto grid = Grid{"O32"};

    IndexKDTree search(geometry());
    search.reserve(grid.size());
    idx_t n{0};
    for (auto& point : grid.lonlat()) {
        search.insert(point, n++);
    }
    // Forgot to call search.build() --> assertion thrown when trying to access
    EXPECT_THROWS_AS(search.closestPoint(PointLonLat{180., 45.}), eckit::AssertionFailed);
}

CASE("test no assertion") {
    // Like case "test assertion", but without reserving size
    auto grid = Grid{"O32"};

    IndexKDTree search(geometry());
    // No search.reserve() --> build() will not be necessary.
    EXPECT(search.empty());
    idx_t n{0};
    for (auto& point : grid.lonlat()) {
        search.insert(point, n++);
        if (ECKIT_515_implemented) {
            EXPECT_EQ(search.size(), n);
        }
    }
    EXPECT_EQ(search.size(), grid.size());
    // search.build() Not required
    EXPECT_NO_THROW(search.closestPoint(PointLonLat{180., 45.}));
}

CASE("test kdtree building with separate lon and lat and payload arrays") {
    IndexKDTree search(geometry());
    search.build(test_lon(), test_lat(), test_payloads());
    validate(search);
}

CASE("test kdtree building with separate lon and lat and raw payload iterators") {
    IndexKDTree search(geometry());
    auto lon       = test_lon();
    auto lat       = test_lat();
    auto payloads_ = test_payloads();
    search.build(lon.begin(), lon.end(), lat.begin(), lat.end(), payloads_.begin(), payloads_.end());
    validate(search);
}

CASE("test kdtree building with separate PointLonLat and payload containers") {
    IndexKDTree search(geometry());
    search.build(test_lonlat(), test_payloads());
    validate(search);
}

CASE("test kdtree building with separate PointXYZ and payload containers") {
    IndexKDTree search(geometry());
    search.build(test_xyz(), test_payloads());
    validate(search);
}

CASE("test assignment") {
    IndexKDTree search;
    search = IndexKDTree();
    search.build(test_lonlat(), test_payloads());
    validate(search);
}

CASE("test closestPoint") {
    auto neighbour          = search().closestPoint(PointLonLat{180., 45.}).payload();
    auto expected_neighbour = 760;
    EXPECT_EQ(neighbour, expected_neighbour);
}

CASE("test closestPoints") {
    auto neighbours          = search().closestPoints(PointLonLat{180., 45.}, 4).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761};
    EXPECT_EQ(neighbours, expected_neighbours);
}

CASE("test closestPointsWithinRadius") {
    double km                = 1000. * radius() / util::Earth::radius();
    auto neighbours          = search().closestPointsWithinRadius(PointLonLat{180., 45.}, 500 * km).payloads();
    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761, 841, 843, 682};
    EXPECT_EQ(neighbours, expected_neighbours);
}

CASE("test compatibility with external eckit KDTree") {
    // External world
    struct ExternalKDTreeTraits {
        using Point   = Point3;
        using Payload = size_t;
    };
    using ExternalKDTree = typename eckit::KDTreeMemory<ExternalKDTreeTraits>;
    auto external_kdtree = std::make_shared<ExternalKDTree>();

    // Atlas world
    IndexKDTree search(external_kdtree, geometry());

    // Construction of tree; could happen in Atlas world or External world separately
    auto grid = Grid{"O32"};
    search.build(grid.lonlat(), PayloadGenerator(grid.size()));

    // Usage in External world
    {
        auto neighbours = external_kdtree->kNearestNeighbours(make_xyz({180., 45.}), 4);
        auto payloads   = std::vector<ExternalKDTree::Payload>{};
        for (auto& neighbour : neighbours) {
            payloads.push_back(neighbour.payload());
        }
        auto expected_payloads = std::vector<ExternalKDTree::Payload>{760, 842, 759, 761};
        EXPECT_EQ(payloads, expected_payloads);
    }

    // Usage in Atlas world
    {
        auto payloads          = search.closestPoints(PointLonLat{180., 45.}, 4).payloads();
        auto expected_payloads = std::vector<IndexKDTree::Payload>{760, 842, 759, 761};
        EXPECT_EQ(payloads, expected_payloads);
    }
}

CASE("test IndexKDTree 2D vs 3D") {
    IndexKDTree2D search2d(geometry());
    IndexKDTree3D search3d(geometry());
    search2d.build(test_lonlat(), test_payloads());
    search3d.build(test_lonlat(), test_payloads());
    auto payloads2d = search2d.closestPoints(PointLonLat{89.9, 44.9}, 4).payloads();
    auto payloads3d = search3d.closestPoints(PointLonLat{89.9, 44.9}, 4).payloads();
    EXPECT_EQ(payloads2d, (std::vector<idx_t>{2, 3, 1, 4}));
    EXPECT_EQ(payloads3d, (std::vector<idx_t>{2, 1, 3, 0}));
    // Note that the expected values are different whether 2D search or 3D search is used
}


CASE("test kdtree with configured geometry") {
    auto grid = Grid{"O32"};

    IndexKDTree search_unit(util::Config("geometry","UnitSphere"));
    IndexKDTree search_earth(util::Config("geometry","Earth"));

    search_unit.build(grid.lonlat(),PayloadGenerator(grid.size()));
    search_earth.build(grid.lonlat(),PayloadGenerator(grid.size()));

    double km_unit                = 1000. / util::Earth::radius();
    double km_earth               = 1000.;
    auto neighbours_unit          = search_unit. closestPointsWithinRadius(PointLonLat{180., 45.}, 500 * km_unit ).payloads();
    auto neighbours_earth         = search_earth.closestPointsWithinRadius(PointLonLat{180., 45.}, 500 * km_earth).payloads();

    auto expected_neighbours = std::vector<idx_t>{760, 842, 759, 761, 841, 843, 682};

    EXPECT_EQ(neighbours_unit,  expected_neighbours);
    EXPECT_EQ(neighbours_earth, expected_neighbours);
}

//------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
