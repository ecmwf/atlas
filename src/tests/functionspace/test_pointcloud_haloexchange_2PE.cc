/*
 * (C) Copyright 2023 ECMWF
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
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


double innerproductwithhalo(const Field& f1, const Field& f2) {
    long sum(0);

    auto view1 = array::make_view<double, 2>(f1);
    auto view2 = array::make_view<double, 2>(f2);

    for (idx_t jn = 0; jn < f1.shape(0); ++jn) {
        for (idx_t jl = 0; jl < f1.levels(); ++jl) {
            sum += view1(jn, jl) * view2(jn, jl);
        }
    }

    atlas::mpi::comm().allReduceInPlace(sum, eckit::mpi::sum());
    return sum;
}


CASE("test_halo_exchange_01") { 
  // Here is some ascii art to describe the test.
  // Remote index is the same for both PEs in this case
  //
  //     _1____0__
  //   1| 0    1 | 0
  //   3| 2    3 | 2
  //    ---------
  //      3    2
  //
  // Point order (locally is also the same)
  //
  //     _4____5__
  //  11| 0    1 | 6
  //  10| 2    3 | 7
  //    ---------
  //      9    8
  //
  // Partition index
  // PE 0:                                       PE 1:
  //
  //     _1____1__                               _0____0__
  //   1| 0    0 | 1                           0| 1    1 | 0
  //   1| 0    0 | 1                           0| 1    1 | 0
  //    ---------                                ---------
  //      1    1                                  0    0
  //
  // Initial Values
  // PE 0:                                       PE 1:
  //
  //     _0____0__                               _0____0__
  //   0|10   11 | 0                           0|20   21 | 0
  //   0|12   13 | 0                           0|22   23 | 0
  //    ---------                                ---------
  //      0    0                                  0    0
  //
  // Values after halo exchange
  // PE 0:                                       PE 1:
  //
  //     21___20__                              11___10__
  //  21|10   11 | 20                        11|20   21 | 10
  //  23|12   13 | 22                        13|22   23 | 12
  //    ---------                               ---------
  //     23   22                                13   12
  //
  // Values after halo exchange and adjoint halo exchange
  // PE 0:                                       PE 1:
  //
  //     _0____0__                              _0___0__
  //   0|30   33 | 0                          0|60   63 | 0
  //   0|36   39 | 0                          0|66   69 | 0
  //    ---------                               ---------
  //      0    0                                 0    0

  double tolerance(1e-16);
  Field lonlat("lonlat", array::make_datatype<double>(), array::make_shape(12, 2));
  Field ghost("ghost", array::make_datatype<int>(), array::make_shape(12));
  Field remote_index("remote_index", array::make_datatype<idx_t>(), array::make_shape(12));
  Field partition("partition", array::make_datatype<int>(), array::make_shape(12));

  auto lonlatv = array::make_view<double, 2>(lonlat);
  auto ghostv = array::make_view<int, 1>(ghost);
  auto remote_indexv = array::make_indexview<idx_t, 1>(remote_index);
  auto partitionv = array::make_view<int, 1>(partition);

  ghostv.assign({ 0, 0, 0, 0,
                  1, 1, 1, 1, 1, 1, 1, 1});

  remote_indexv.assign({0, 1, 2, 3, 1, 0, 0, 2, 2, 3, 3, 1});


  if (atlas::mpi::rank() == 0) {
    // center followed by clockwise halo starting from top left
    lonlatv.assign({-45.0,  45.0,  45.0,  45.0,    // center, first row
                    -45.0, -45.0,  45.0, -45.0,    // center, second row
                    225.0,  45.0, 135.0,  45.0,    // up
                    135.0,  45.0, 135.0, -45.0,    // right
                    135.0, -45.0, 225.0, -45.0,    // down
                    225.0, -45.0, 225.0,  45.0});  // left

    partitionv.assign({0, 0, 0, 0,
                       1, 1, 1, 1, 1, 1, 1, 1});



  } else if (atlas::mpi::rank() == 1) {
    // center followed by clockwise halo starting from top left
    lonlatv.assign({135.0,  45.0, 225.0,  45.0,
                    135.0, -45.0, 225.0, -45.0,
                     45.0,  45.0, -45.0,  45.0,
                    -45.0,  45.0, -45.0, -45.0,
                    -45.0, -45.0,  45.0, -45.0,
                     45.0, -45.0,  45.0,  45.0});

    partitionv.assign({1, 1, 1, 1,
                       0, 0, 0, 0, 0, 0, 0, 0});

  }

  FieldSet fset;
  fset.add(lonlat);
  fset.add(ghost);
  fset.add(remote_index);
  fset.add(partition);

  auto fs2 = functionspace::PointCloud(fset);
  Field f1 = fs2.createField<double>(option::name("f1") | option::levels(2));
  Field f2 = fs2.createField<double>(option::name("f2") | option::levels(2));
  auto f1v = array::make_view<double, 2>(f1);
  auto f2v = array::make_view<double, 2>(f2);

  f1v.assign(0.0);
  f2v.assign(0.0);
  for (idx_t i = 0; i < f2v.shape(0); ++i) {
    for (idx_t l = 0; l < f2v.shape(1); ++l) {
      auto ghostv2 = array::make_view<int, 1>(fs2.ghost());
      if (ghostv2(i) == 0) {
        f1v(i, l) = (atlas::mpi::rank() +1) * 10.0  + i;
        f2v(i, l) =  f1v(i, l);
      }
    }
  }

  f2.haloExchange();

  // adjoint test
  double sum1 = innerproductwithhalo(f2, f2);

  f2.adjointHaloExchange();

  double sum2 = innerproductwithhalo(f1, f2);

  atlas::mpi::comm().allReduceInPlace(sum1, eckit::mpi::sum());
  atlas::mpi::comm().allReduceInPlace(sum2, eckit::mpi::sum());
  EXPECT(std::abs(sum1 - sum2)/ std::abs(sum1) < tolerance);
  atlas::Log::info() << "adjoint test passed :: "
                     << "sum1 " << sum1 << " sum2 " << sum2 << " normalised difference "
                     << std::abs(sum1 - sum2)/ std::abs(sum1) << std::endl;

  // In this case the effect of the halo exchange followed by the adjoint halo exchange
  // multiples the values by a factor of 3 (see pictures above)
  for (idx_t i = 0; i < f2v.shape(0); ++i) {
    for (idx_t l = 0; l < f2v.shape(1); ++l) {
       EXPECT( std::abs(f2v(i, l) - 3.0 * f1v(i, l)) < tolerance);
    }
  }
  atlas::Log::info() << "values from halo followed by halo adjoint are as expected "
                     << std::endl;
}


CASE("test_halo_exchange_02") { 

  double tolerance(1e-16);

  Field lonlat("lonlat", array::make_datatype<double>(), array::make_shape(12, 2));
  Field ghost("ghost", array::make_datatype<int>(), array::make_shape(12));

  auto lonlatv = array::make_view<double, 2>(lonlat);
  auto ghostv = array::make_view<int, 1>(ghost);

  if (atlas::mpi::rank() == 0) {
    // center followed by clockwise halo starting from top left
    lonlatv.assign({-45.0,  45.0,  45.0,  45.0,    // center, first row
                    -45.0, -45.0,  45.0, -45.0,    // center, second row
                    225.0,  45.0, 135.0,  45.0,    // up
                    135.0,  45.0, 135.0, -45.0,    // right
                    135.0, -45.0, 225.0, -45.0,    // down
                    225.0, -45.0, 225.0,  45.0});  // left

    ghostv.assign({ 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1});

  } else if (atlas::mpi::rank() == 1) {
    // center followed by clockwise halo starting from top left
    lonlatv.assign({135.0,  45.0, 225.0,  45.0,
                    135.0, -45.0, 225.0, -45.0,
                     45.0,  45.0, -45.0,  45.0,
                    -45.0,  45.0, -45.0, -45.0,
                    -45.0, -45.0,  45.0, -45.0,
                     45.0, -45.0,  45.0,  45.0});

    ghostv.assign({ 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1});
  }


  auto pcfs = functionspace::PointCloud(lonlat, ghost);

  Field f1 = pcfs.createField<double>(option::name("f1") | option::levels(2));
  Field f2 = pcfs.createField<double>(option::name("f2") | option::levels(2));
  auto f1v = array::make_view<double, 2>(f1);
  auto f2v = array::make_view<double, 2>(f2);

  f1v.assign(0.0);
  f2v.assign(0.0);

  for (idx_t i = 0; i < f2v.shape(0); ++i) {
    for (idx_t l = 0; l < f2v.shape(1); ++l) {
      auto ghostv2 = array::make_view<int, 1>(pcfs.ghost());
      if (ghostv2(i) == 0) {
        f1v(i, l) = (atlas::mpi::rank() +1) * 10.0  + i;
        f2v(i, l) =  f1v(i, l);
      }
    }
  }

  f2.haloExchange();

   // adjoint test
  double sum1 = innerproductwithhalo(f2, f2);

  f2.adjointHaloExchange();

  double sum2 = innerproductwithhalo(f1, f2);

  atlas::mpi::comm().allReduceInPlace(sum1, eckit::mpi::sum());
  atlas::mpi::comm().allReduceInPlace(sum2, eckit::mpi::sum());
  EXPECT(std::abs(sum1 - sum2)/ std::abs(sum1) < tolerance);
  atlas::Log::info() << "adjoint test passed :: "
                     << "sum1 " << sum1 << " sum2 " << sum2 << " normalised difference "
                     << std::abs(sum1 - sum2)/ std::abs(sum1) << std::endl;

  // In this case the effect of the halo exchange followed by the adjoint halo exchange
  // multiples the values by a factor of 3 (see pictures above)
  for (idx_t i = 0; i < f2v.shape(0); ++i) {
    for (idx_t l = 0; l < f2v.shape(1); ++l) {
       EXPECT( std::abs(f2v(i, l) - 3.0 * f1v(i, l)) < tolerance);
    }
  }
  atlas::Log::info() << "values from halo followed by halo adjoint are as expected "
                     << std::endl;

}


}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
