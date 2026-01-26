/*
 * (C) Copyright 2023 ECMWF
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

using namespace eckit;
using namespace atlas::functionspace;
using namespace atlas::util;

namespace atlas {
namespace test {


CASE("test_halo_exchange_01") {

  //
  //     ++ point order ++
  //
  //       PE 0     PE 1
  //     _________________
  //    | 0    1 || 0   1 |
  //    | 3    2 || 3   2 |
  //    ------------------
  //    ------------------
  //    | 0    1    2   3 |    
  //    | 7    6    5   4 |
  //    ------------------
  //           PE 2
  //
  // 'own points': clockwise starting from top left
  //


  // number of points (own + ghost)
  int no_points{0};
  if (atlas::mpi::rank() == 0) {
    no_points = 8;
  } else if (atlas::mpi::rank() == 1) {
    no_points = 8;
  } else if (atlas::mpi::rank() == 2) {
    no_points = 12;
  }

  Field lonlat("lonlat", array::make_datatype<double>(),
    array::make_shape(no_points, 2));
  // ghost points flags: 0={is not a ghost point}, 1={is a ghost point}
  Field gpoints("ghost", array::make_datatype<int>(),
    array::make_shape(no_points));

  auto lonlatv = array::make_view<double, 2>(lonlat);
  auto gpointsv = array::make_view<int, 1>(gpoints);

  if (atlas::mpi::rank() == 0) {
    // own points: clockwise starting from top left
    // halo points: clockwise starting from top left
    lonlatv.assign({0.0,  1.0,  1.0,  1.0,    // center, first row  [own]
                    1.0, -1.0,  0.0, -1.0,    // center, second row [own]
              	    2.0,  1.0,  2.0, -1.0,    // right              [ghost]
                    1.0, -2.0,  0.0, -2.0});  // down               [ghost]
 
    gpointsv.assign({0, 0, 0, 0,
                     1, 1, 1, 1});

  } else if (atlas::mpi::rank() == 1) {
    lonlatv.assign({2.0,  1.0,  3.0,  1.0,    // center, first row  [own]
                    3.0, -1.0,  2.0, -1.0,    // center, second row [own]
                    3.0, -2.0,  2.0, -2.0,    // down               [ghost]
                    1.0, -1.0,  1.0,  1.0});  // left               [ghost]
 
    gpointsv.assign({0, 0, 0, 0,
                     1, 1, 1, 1});

  } else if (atlas::mpi::rank() == 2) {
    lonlatv.assign({0.0, -2.0, 1.0, -2.0, 2.0, -2.0, 3.0, -2.0,   // center, first row  [own]
                    3.0, -3.0, 2.0, -3.0, 1.0, -3.0, 0.0, -3.0,   // center, second row [own]
                    0.0, -1.0, 1.0, -1.0, 2.0, -1.0, 3.0, -1.0}); // top                [ghost]
 
    gpointsv.assign({0, 0, 0, 0, 0, 0, 0, 0,
                     1, 1, 1, 1});

  }


  // function space
  auto pcfs = functionspace::PointCloud(lonlat, gpoints);

  // remote indexes (reference)
  std::vector<idx_t> remote_idxs_ref;
  if (atlas::mpi::rank() == 0) {
    remote_idxs_ref = {0, 1, 2, 3, 0, 3, 1, 0};
  } else if (atlas::mpi::rank() == 1) {
    remote_idxs_ref = {0, 1, 2, 3, 3, 2, 2, 1};
  } else if (atlas::mpi::rank() == 2) {
    remote_idxs_ref = {0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 3, 2};
  }

  // remote indexes
  auto remote_index = pcfs.remote_index();

  EXPECT(remote_index.size() == (lonlat.size()/2));;

  auto remote_indexv = array::make_indexview<idx_t, 1>(remote_index);
  for (idx_t i = 0; i < remote_indexv.shape(0); ++i) {
    EXPECT(remote_indexv(i) == remote_idxs_ref.at(i));;
  }

}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
