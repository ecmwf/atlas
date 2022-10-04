/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"

#include "tests/AtlasTestEnvironment.h"

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------

std::vector<PointXY> ref_lonlat() {
    static std::vector<PointXY> v = {{00., 0.}, {10., 0.}, {20., 0.}, {30., 0.}, {40., 0.},
                                     {50., 0.}, {60., 0.}, {70., 0.}, {80., 0.}, {90., 0.}};
    return v;
}
std::vector<PointXY> ref_xy() {
    static std::vector<PointXY> v = {{00., 0.}, {10., 0.}, {20., 0.}, {30., 0.}, {40., 0.},
                                     {50., 0.}, {60., 0.}, {70., 0.}, {80., 0.}, {90., 0.}};
    return v;
}
std::vector<PointXYZ> ref_xyz() {
    static std::vector<PointXYZ> v = {{00., 0., 1.}, {10., 0., 2.}, {20., 0., 3.}, {30., 0., 4.}, {40., 0., 5.},
                                      {50., 0., 6.}, {60., 0., 7.}, {70., 0., 8.}, {80., 0., 9.}, {90., 0., 10.}};
    return v;
}

template <typename C1, typename C2>
bool equal(const C1& c1, const C2& c2) {
    return c1.size() == c2.size() && std::equal(std::begin(c1), std::end(c1), std::begin(c2));
}


CASE("test_functionspace_PointCloud create from field") {
    Field points("points", array::make_datatype<double>(), array::make_shape(10, 2));
    auto xy = array::make_view<double, 2>(points);
    xy.assign({00., 0., 10., 0., 20., 0., 30., 0., 40., 0., 50., 0., 60., 0., 70., 0., 80., 0., 90., 0.});

    functionspace::PointCloud pointcloud(points);
    EXPECT(pointcloud.size() == 10);

    points.dump(Log::info());
    Log::info() << std::endl;

    EXPECT_THROWS(pointcloud.iterate().xyz());
    EXPECT(equal(pointcloud.iterate().xy(), ref_xy()));
    EXPECT(equal(pointcloud.iterate().lonlat(), ref_lonlat()));
}

CASE("test_functionspace_PointCloud create from 2D vector") {
    std::vector<PointXY> points{{00., 0.}, {10., 0.}, {20., 0.}, {30., 0.}, {40., 0.},
                                {50., 0.}, {60., 0.}, {70., 0.}, {80., 0.}, {90., 0.}};

    functionspace::PointCloud pointcloud(points);
    EXPECT(pointcloud.size() == 10);

    pointcloud.lonlat().dump(Log::info());
    Log::info() << std::endl;

    EXPECT_THROWS(pointcloud.iterate().xyz());
    EXPECT(equal(pointcloud.iterate().xy(), ref_xy()));
    EXPECT(equal(pointcloud.iterate().lonlat(), ref_lonlat()));
}

CASE("test_functionspace_PointCloud create from 3D vector") {
    std::vector<PointXYZ> points{{00., 0., 1.}, {10., 0., 2.}, {20., 0., 3.}, {30., 0., 4.}, {40., 0., 5.},
                                 {50., 0., 6.}, {60., 0., 7.}, {70., 0., 8.}, {80., 0., 9.}, {90., 0., 10.}};
    functionspace::PointCloud pointcloud(points);
    EXPECT(pointcloud.size() == 10);

    pointcloud.lonlat().dump(Log::info());
    Log::info() << std::endl;

    EXPECT(equal(pointcloud.iterate().xyz(), ref_xyz()));
    EXPECT(equal(pointcloud.iterate().xy(), ref_xy()));
    EXPECT(equal(pointcloud.iterate().lonlat(), ref_lonlat()));
}


CASE("test_functionspace_PointCloud create from 2D initializer list") {
    functionspace::PointCloud pointcloud = {{00., 0.}, {10., 0.}, {20., 0.}, {30., 0.}, {40., 0.},
                                            {50., 0.}, {60., 0.}, {70., 0.}, {80., 0.}, {90., 0.}};
    EXPECT(pointcloud.size() == 10);

    pointcloud.lonlat().dump(Log::info());
    Log::info() << std::endl;

    EXPECT_THROWS(pointcloud.iterate().xyz());
    EXPECT(equal(pointcloud.iterate().xy(), ref_xy()));
    EXPECT(equal(pointcloud.iterate().lonlat(), ref_lonlat()));
}

CASE("test_functionspace_PointCloud create from 3D initializer list") {
    functionspace::PointCloud pointcloud = {{00., 0., 1.}, {10., 0., 2.}, {20., 0., 3.}, {30., 0., 4.}, {40., 0., 5.},
                                            {50., 0., 6.}, {60., 0., 7.}, {70., 0., 8.}, {80., 0., 9.}, {90., 0., 10.}};
    EXPECT(pointcloud.size() == 10);

    pointcloud.lonlat().dump(Log::info());
    Log::info() << std::endl;
    pointcloud.vertical().dump(Log::info());
    Log::info() << std::endl;

    EXPECT(equal(pointcloud.iterate().xyz(), ref_xyz()));
    EXPECT(equal(pointcloud.iterate().xy(), ref_xy()));
    EXPECT(equal(pointcloud.iterate().lonlat(), ref_lonlat()));
}


//-----------------------------------------------------------------------------

CASE("test_createField") {
    FunctionSpace p1;
    {
        Field points("points", array::make_datatype<double>(), array::make_shape(10, 2));
        auto xy = array::make_view<double, 2>(points);
        xy.assign({00., 0., 10., 0., 20., 0., 30., 0., 40., 0., 50., 0., 60., 0., 70., 0., 80., 0., 90., 0.});
        p1 = functionspace::PointCloud(points);
    }

    Field f1 = p1.createField<double>(option::name("f1") | option::levels(3));
    EXPECT_EQ(f1.levels(), 3);
    EXPECT_EQ(f1.shape(0), 10);
    EXPECT_EQ(f1.shape(1), 3);

    FunctionSpace p2;
    {
        Field points("points", array::make_datatype<double>(), array::make_shape(4, 2));
        auto xy = array::make_view<double, 2>(points);
        xy.assign({20., 0., 40., 0., 70., 0., 90., 0.});
        p2 = functionspace::PointCloud(points);
    }
    Field f2 = p2.createField(f1, util::NoConfig());
    EXPECT_EQ(f2.levels(), 3);
    EXPECT_EQ(f2.shape(0), 4);
    EXPECT_EQ(f2.shape(1), 3);

    Field f3 = p2.createField(f1, option::variables(5));
    EXPECT_EQ(f3.levels(), 3);
    EXPECT_EQ(f3.shape(0), 4);
    EXPECT_EQ(f3.shape(1), 3);
    EXPECT_EQ(f3.shape(2), 5);
}

CASE("test_createFieldSet") {

  Field lonlat("lonlat", array::make_datatype<double>(), array::make_shape(12, 2));
  Field ghost("ghost", array::make_datatype<int>(), array::make_shape(12));
  Field remote_index("remote_index", array::make_datatype<int>(), array::make_shape(12));
  Field partition("partition", array::make_datatype<int>(), array::make_shape(12));

  auto lonlatv = array::make_view<double, 2>(lonlat);
  auto ghostv = array::make_view<int, 1>(ghost);
  auto remote_indexv = array::make_view<int, 1>(remote_index);
  auto partitionv = array::make_view<int, 1>(partition);

  ghostv.assign({ 0, 0, 0, 0,
                  1, 1, 1, 1, 1, 1, 1, 1});

  remote_indexv.assign({0, 1, 2, 3, 1, 0, 0, 2, 2, 3, 1});


  if (atlas::mpi::rank() == 0) {
    // center followed by clockwise halo starting from top left
    lonlatv.assign({-45.0, 45.0,
                    45.0,  45.0,
                   -45.0, -45.0,
                    -45.0, 45.0,  // center
                    225.0, 45.0, 135.0, 45.0,  // up
                    135.0, 45.0, 135.0, -45.0,  //left
                    135.0, -45.0, 135.0, 45.0,  // down
                    135.0, 45.0,   225.0, 45.0}); // right

    partitionv.assign({0, 0, 0, 0,
                       1, 1, 1, 1, 1, 1, 1, 1});



  } else if (atlas::mpi::rank() == 1) {
    // center followed by clockwise halo starting from top left
    lonlatv.assign({135.0, 45.0, 225.0, 45.0, 135.0, -45.0, 135.0, 45.0, // center
                    45.0, 45.0, -45.0, 45.0, // up
                   -45.0, 45.0, -45.0, -45.0, // left
                   -45.0, -45.0, -45.0, 45.0, // down
                   -45.0, 45.0,   45.0, 45.0}); // right

    partitionv.assign({1, 1, 1, 1,
                       0, 0, 0, 0, 0, 0, 0, 0});

  }

  atlas::FieldSet fset;
  fset.add(lonlat);
  fset.add(ghost);
  fset.add(remote_index);
  fset.add(partition);

  auto fs2 = functionspace::PointCloud(fset);
  Field f2 = fs2.createField<double>(option::name("f1") | option::levels(2) | option::halo(1));
  auto f2v = array::make_view<double, 2>(f2);

  f2v.assign(0.0);
  for (idx_t i = 0; i < f2v.shape(0); ++i) {
    for (idx_t l = 0; l < f2v.shape(1); ++l) {
      auto ghostv2 = array::make_view<int, 1>(fs2.ghost());
      if (ghostv2(i) == 0) {
        f2v(i, l) = (atlas::mpi::rank() +1) * 10.0  + i;
      }
    }
  }

  f2.haloExchange();

  for (idx_t i = 0; i < f2v.shape(0); ++i) {
    for (idx_t l = 0; l < f2v.shape(1); ++l) {
      std::cout << "i l rank value = " << i << " " << l << " " << atlas::mpi::rank() << " " << f2v(i, l) << std::endl;
    }
  }



}



//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
