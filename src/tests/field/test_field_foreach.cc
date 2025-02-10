/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "tests/AtlasTestEnvironment.h"

#include "atlas/field/for_each.h"
#include "atlas/option.h"

using namespace atlas::array;
using namespace atlas::array::helpers;

namespace atlas {
namespace test {

auto split_index_3d = [] (int idx, auto shape) {
  int i = idx / (shape[2]*shape[1]);
  int j = (idx - i * (shape[2] * shape[1]))/shape[2];
  int k = idx - i * (shape[2]*shape[1]) - j*shape[2];
  return std::make_tuple(i,j,k);
};

CASE( "test field::for_each_value" ) {
  Field f("name", array::make_datatype<double>(), array::make_shape({4,2,8}));

  int v=0;
  field::for_each_value(f, [&](double& x) {
    auto [i,j,k] = split_index_3d(v,f.shape());
    x = 100*i + 10*j + k;
    v++;
  });

  v = 0;
  field::for_each_value(f, [&](double& x) {
    auto [i,j,k] = split_index_3d(v,f.shape());
    double expected = 100*i + 10*j + k;
    EXPECT_EQ(x, expected);
    ++v;
  });
}


CASE( "test field::for_each_value_masked; horizontal_dimension {0} (implicit)" ) {
  Field f("name", array::make_datatype<double>(), array::make_shape({4,2,8}));

  array::make_view<double,3>(f).assign(0.);

  Field ghost("ghost", array::make_datatype<int>(), array::make_shape(4) );
  array::make_view<int,1>(ghost).assign(0);
  array::make_view<int,1>(ghost)(2) = 1;

  field::for_each_value_masked(ghost, f, [&](double& x) {
    x = 1.;
  });

  int v = 0;
  field::for_each_value(f, [&](double& x) {
    auto [i,j,k] = split_index_3d(v,f.shape());
    double expected = 1.;
    if (i==2) {
      expected = 0.;
    }
    EXPECT_EQ(x, expected);
    ++v;
  });
}

CASE( "test field::for_each_value_masked; horizontal_dimension {0,2} mask_rank1" ) {
  Field f("name", array::make_datatype<double>(), array::make_shape({4,2,8}));
  f.set_horizontal_dimension({0,2});

  array::make_view<double,3>(f).assign(0.);

  Field ghost("ghost", array::make_datatype<int>(), array::make_shape(f.shape(0)*f.shape(2)) );
  array::make_view<int,1>(ghost).assign(0);
  array::make_view<int,1>(ghost)(2*f.shape(2)+3) = 1;

  field::for_each_value_masked(ghost, f, [&](double& x) {
    x = 1.;
  });

  int v = 0;
  field::for_each_value(f, [&](double& x) {
    auto [i,j,k] = split_index_3d(v,f.shape());
    double expected = 1.;
    if (i==2 && k==3) {
      expected = 0.;
    }
    EXPECT_EQ(x, expected);
    ++v;
  });
}

CASE( "test field::for_each_value_masked;  horizontal_dimension {0,2} mask_rank2" ) {
  Field f("name", array::make_datatype<double>(), array::make_shape({4,2,8}));
  f.set_horizontal_dimension({0,2});

  array::make_view<double,3>(f).assign(0.);

  Field ghost("ghost", array::make_datatype<int>(), array::make_shape(f.shape(0),f.shape(2)) );
  array::make_view<int,2>(ghost).assign(0);
  array::make_view<int,2>(ghost)(2,3) = 1;

  field::for_each_value_masked(ghost, f, [&](double& x) {
    x = 1.;
  });

  int v = 0;
  field::for_each_value(f, [&](double& x) {
    auto [i,j,k] = split_index_3d(v,f.shape());
    double expected = 1.;
    if (i==2 && k==3) {
      expected = 0.;
    }
    EXPECT_EQ(x, expected);
    ++v;
  });
}

CASE( "test field::for_each_column" ) {
  Field f("name", array::make_datatype<double>(), array::make_shape({6,4,2}));

  int v=0;
  field::for_each_value(execution::seq, f, [&](double& x) {
    auto [i,j,k] = split_index_3d(v,f.shape());
    x = 100*i + 10*j + k;
    v++;
  });

  [[maybe_unused]] auto print_column_1d =  [&](const array::View<double,1>& column) {
      Log::info() << "";
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        if (jlev != 0 ) {
          Log::info() << ", ";
        }
        Log::info() << column(jlev);
      }
      Log::info() << "," << std::endl;
  };

  [[maybe_unused]] auto print_column_2d = [&](array::View<double,2>& column) {
    Log::info() << "";
    for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
      if (jlev != 0 ) {
        Log::info() << ",\n";
      }
      for( idx_t jvar=0; jvar<column.shape(1); ++jvar ) {
        if (jvar != 0 ) {
          Log::info() << ", ";
        }
        Log::info() << column(jlev,jvar);
      }
    }
    Log::info() << "," << std::endl;
  };

  SECTION( "horizontal_dimension = {0,1}") {
    f.set_horizontal_dimension({0,1});

    std::vector<double> visited_values;
    field::for_each_column(f, [&visited_values](const array::View<double,1>& column) {
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        visited_values.push_back(column(jlev));
      }
    });
    EXPECT_EQ( visited_values.size(), f.size() );
    std::vector<double> expected_values = {
      0, 1,     //
      10, 11,   //
      20, 21,   //
      30, 31,   //
      100, 101, //
      110, 111, //
      120, 121, //
      130, 131, //
      200, 201, //
      210, 211, //
      220, 221, //
      230, 231, //
      300, 301, //
      310, 311, //
      320, 321, //
      330, 331, //
      400, 401, //
      410, 411, //
      420, 421, //
      430, 431, //
      500, 501, //
      510, 511, //
      520, 521, //
      530, 531  //
    };
    EXPECT_EQ( visited_values, expected_values );
  }
  SECTION( "horizontal_dimension = {0,2}") {
    f.set_horizontal_dimension({0,2});

    std::vector<double> visited_values;
    field::for_each_column(f, [&visited_values](const array::View<double,1>& column) {
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        visited_values.push_back(column(jlev));
      }
    });
    EXPECT_EQ( visited_values.size(), f.size() );
    std::vector<double> expected_values = {
      0, 10, 20, 30,      //
      1, 11, 21, 31,      //
      100, 110, 120, 130, //
      101, 111, 121, 131, //
      200, 210, 220, 230, //
      201, 211, 221, 231, //
      300, 310, 320, 330, //
      301, 311, 321, 331, //
      400, 410, 420, 430, //
      401, 411, 421, 431, // 
      500, 510, 520, 530, //
      501, 511, 521, 531  //
    };
    EXPECT_EQ( visited_values, expected_values );
  }
  SECTION( "horizontal_dimension = {1,2}") {
    f.set_horizontal_dimension({1,2});

    std::vector<double> visited_values;
    field::for_each_column(f, [&visited_values](const array::View<double,1>& column) {
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        visited_values.push_back(column(jlev));
      }
    });
    // field::for_each_column(execution::seq, f, print_column_1d);
    EXPECT_EQ( visited_values.size(), f.size() );
    std::vector<double> expected_values = {
      0, 100, 200, 300, 400, 500,  //
      1, 101, 201, 301, 401, 501,  //
      10, 110, 210, 310, 410, 510, //
      11, 111, 211, 311, 411, 511, //
      20, 120, 220, 320, 420, 520, //
      21, 121, 221, 321, 421, 521, //
      30, 130, 230, 330, 430, 530, //
      31, 131, 231, 331, 431, 531  //
    };
    EXPECT_EQ( visited_values, expected_values );
  }

  SECTION( "horizontal_dimension = {0}") {
    f.set_horizontal_dimension({0});

    std::vector<double> visited_values;
    field::for_each_column(f, [&visited_values](const array::View<double,2>& column) {
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        for( idx_t jvar=0; jvar<column.shape(1); ++jvar) {
          visited_values.push_back(column(jlev,jvar));
        }
      }
    });
    EXPECT_EQ( visited_values.size(), f.size() );
        // field::for_each_column(execution::seq, f, print_column_2d);

    std::vector<double> expected_values = {
      0, 1,     //
      10, 11,   //
      20, 21,   //
      30, 31,   //
      100, 101, //
      110, 111, //
      120, 121, //
      130, 131, //
      200, 201, //
      210, 211, //
      220, 221, //
      230, 231, //
      300, 301, //
      310, 311, //
      320, 321, //
      330, 331, //
      400, 401, //
      410, 411, //
      420, 421, //
      430, 431, //
      500, 501, //
      510, 511, //
      520, 521, //
      530, 531  //
    };
    EXPECT_EQ( visited_values, expected_values );
  }

  SECTION( "horizontal_dimension = {1}") {
    f.set_horizontal_dimension({1});

    std::vector<double> visited_values;
    field::for_each_column(f, [&visited_values](const array::View<double,2>& column) {
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        for( idx_t jvar=0; jvar<column.shape(1); ++jvar) {
          visited_values.push_back(column(jlev,jvar));
        }
      }
    });
    EXPECT_EQ( visited_values.size(), f.size() );

    std::vector<double> expected_values = {
      0, 1,     //
      100, 101, //
      200, 201, //
      300, 301, //
      400, 401, //
      500, 501, //
      10, 11,   //
      110, 111, //
      210, 211, //
      310, 311, //
      410, 411, //
      510, 511, //
      20, 21,   //
      120, 121, //
      220, 221, //
      320, 321, //
      420, 421, //
      520, 521, //
      30, 31,   //
      130, 131, //
      230, 231, //
      330, 331, //
      430, 431, //
      530, 531  //
    };
    EXPECT_EQ( visited_values, expected_values );
  }

  SECTION( "horizontal_dimension = {2}") {
    f.set_horizontal_dimension({2});

    std::vector<double> visited_values;
    field::for_each_column(f, [&visited_values](const array::View<double,2>& column) {
      for( idx_t jlev=0; jlev<column.shape(0); ++jlev) {
        for( idx_t jvar=0; jvar<column.shape(1); ++jvar) {
          visited_values.push_back(column(jlev,jvar));
        }
      }
    });
    EXPECT_EQ( visited_values.size(), f.size() );
    std::vector<double> expected_values = {
      0, 10, 20, 30,      //
      100, 110, 120, 130, //
      200, 210, 220, 230, //
      300, 310, 320, 330, //
      400, 410, 420, 430, //
      500, 510, 520, 530, //
      1, 11, 21, 31,      //
      101, 111, 121, 131, //
      201, 211, 221, 231, //
      301, 311, 321, 331, //
      401, 411, 421, 431, //
      501, 511, 521, 531  //
    };
    EXPECT_EQ( visited_values, expected_values );
  }
}

CASE( "test field::for_each_column multiple" ) {
  Field f1("f1", array::make_datatype<double>(),array::make_shape({6,4,2}));
  Field f2("f2", array::make_datatype<double>(),array::make_shape({6,4,2}));
  Field f3("f2", array::make_datatype<double>(),array::make_shape({6,4,2}));

  int v=0;
  field::for_each_value(f1, f2, [&](double& x1, double& x2) {
    auto [i,j,k] = split_index_3d(v,f1.shape());
    x1 = 100*i + 10*j + k;
    x2 = 2*x1;
    v++;
  });

  auto check_result = [&]() {
    auto view = array::make_view<double,3>(f3);
    for (idx_t i=0; i<view.shape(0); ++i) {
      for (idx_t j=0; j<view.shape(1); ++j) {
        for (idx_t k=0; k<view.shape(2); ++k) {
          double expected = (100*i + 10*j + k) * 3;
          EXPECT_EQ(view(i,j,k), expected);
        }
      }
    }
  };
  constexpr int any=-1;
  auto check_result_masked = [&](int _i, int _j, int _k) {
    auto view = array::make_view<double,3>(f3);
    for (idx_t i=0; i<view.shape(0); ++i) {
      for (idx_t j=0; j<view.shape(1); ++j) {
        for (idx_t k=0; k<view.shape(2); ++k) {
          double expected = (100*i + 10*j + k) * 3;;
          if ( (i == _i || _i == any) && (j == _j || _j == any) && (k == _k || _k == any) ) {
            expected = 0;
          }
          EXPECT_EQ(view(i,j,k), expected);
        }
      }
    }
  };

  auto make_ghost = [](Field f1, std::vector<int> mask) {
    auto h_dim = f1.horizontal_dimension();
    ATLAS_ASSERT(h_dim.size() == mask.size());
    std::vector<int> h_shape;
    for (auto h: h_dim) {
      h_shape.emplace_back(f1.shape(h));
    }

    size_t h_size = 1;
    for (auto h: h_shape) {
      h_size *= h;
    }

    Field ghost("ghost", array::make_datatype<int>(), array::make_shape(h_size));
    auto ghost_v = array::make_view<int,1>(ghost);
    ghost_v.assign(0);
    if( mask.size() == 1 ) {
      ghost_v(mask[0]) = 1;
    }
    if( mask.size() == 2 ) {
      ghost_v(h_shape[1]*mask[0] + mask[1]) = 1;
    }
    if( mask.size() == 3 ) {
      ghost_v(h_shape[2]*h_shape[1]*mask[0] + h_shape[1]*mask[1] + mask[2]) = 1;
    }
    return ghost;
  };

  SECTION( "for_each_column; horizontal_dimension = {0,1}") {
    f1.set_horizontal_dimension({0,1});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column(f1, f2, f3, [](
        array::View<const double,1>&& c1,
        array::View<const double,1>&& c2,
        array::View<double,1>&& c3) {
      for( idx_t jlev=0; jlev<c1.shape(0); ++jlev) {
        c3(jlev) = c1(jlev) + c2(jlev);
      }
    });
    check_result();
  }

  SECTION( "for_each_column; horizontal_dimension = {0,2}") {
    f1.set_horizontal_dimension({0,2});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column(f1, f2, f3, [](
        array::View<const double,1>&& c1,
        array::View<const double,1>&& c2,
        array::View<double,1>&& c3) {
      for( idx_t jlev=0; jlev<c1.shape(0); ++jlev) {
        c3(jlev) = c1(jlev) + c2(jlev);
      }
    });
    check_result();
  }

  SECTION( "for_each_column; horizontal_dimension = {1,2}") {
    f1.set_horizontal_dimension({1,2});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column(f1, f2, f3, [](
        array::View<const double,1>&& c1,
        array::View<const double,1>&& c2,
        array::View<double,1>&& c3) {
      for( idx_t jlev=0; jlev<c1.shape(0); ++jlev) {
        c3(jlev) = c1(jlev) + c2(jlev);
      }
    });
    check_result();
  }

  SECTION( "for_each_column; horizontal_dimension = {0}") {
    f1.set_horizontal_dimension({0});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column(f1, f2, f3, [](
        array::View<const double,2>&& c1,
        array::View<const double,2>&& c2,
        array::View<double,2>&& c3) {
      for( idx_t jlev=0; jlev<c1.shape(0); ++jlev) {
        for( idx_t jvar=0; jvar<c1.shape(1); ++jvar) {
          c3(jlev,jvar) = c1(jlev,jvar) + c2(jlev,jvar);
        }
      }
    });
    check_result();
  }

  SECTION( "for_each_column; horizontal_dimension = {1}") {
    f1.set_horizontal_dimension({1});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column(f1, f2, f3, [](
        array::View<const double,2>&& c1,
        array::View<const double,2>&& c2,
        array::View<double,2>&& c3) {
      for( idx_t jlev=0; jlev<c1.shape(0); ++jlev) {
        for( idx_t jvar=0; jvar<c1.shape(1); ++jvar) {
          c3(jlev,jvar) = c1(jlev,jvar) + c2(jlev,jvar);
        }
      }
    });
    check_result();
  }

  SECTION( "for_each_column; horizontal_dimension = {2}") {
    f1.set_horizontal_dimension({2});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column(f1, f2, f3, [](
        array::View<const double,2>&& c1,
        array::View<const double,2>&& c2,
        array::View<double,2>&& c3) {
      for( idx_t jlev=0; jlev<c1.shape(0); ++jlev) {
        for( idx_t jvar=0; jvar<c1.shape(1); ++jvar) {
          c3(jlev,jvar) = c1(jlev,jvar) + c2(jlev,jvar);
        }
      }
    });
    check_result();
  }

  SECTION( "for_each_column_masked; horizontal_dimension = {0}") {
    f1.set_horizontal_dimension({0});
    auto ghost = make_ghost(f1, {4});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column_masked(ghost,
      f1, f2, f3, [](
        array::View<const double,2>&& c1,
        array::View<const double,2>&& c2,
        array::View<double,2>&& c3) {
        for( idx_t jlev=0; jlev<c3.shape(0); ++jlev) {
          for( idx_t jvar=0; jvar<c3.shape(1); ++jvar) {
            c3(jlev,jvar) = c1(jlev,jvar) + c2(jlev,jvar);
          }
        }
    });
    check_result_masked(4,any,any);
  }

  SECTION( "for_each_column_masked; horizontal_dimension = {0,2}") {
    f1.set_horizontal_dimension({0,2});
    auto ghost = make_ghost(f1, {3,1});
    array::make_view<double,3>(f3).assign(0.);
    field::for_each_column_masked(ghost, f1, f2, f3, [](
        array::View<const double,1>&& c1,
        array::View<const double,1>&& c2,
        array::View<double,1>&& c3) {
        for( idx_t jlev=0; jlev<c3.shape(0); ++jlev) {
          c3(jlev) = c1(jlev) + c2(jlev);
        }
    });
    check_result_masked(3,any,1);
  }


}


CASE( "test execution policy" ) {
  idx_t ni = 1000;
  idx_t nj = 100;
  idx_t nk = 16;
  // ni /= 4;
  nj /= 4;
  nk /= 4;
  Field f1("f1", array::make_datatype<double>(),array::make_shape({ni,nj,nk}));
  Field f2("f2", array::make_datatype<double>(),array::make_shape({ni,nj,nk}));
  Field f3("f2", array::make_datatype<double>(),array::make_shape({ni,nj,nk}));

  int v=0;
  field::for_each_value(f1, f2, [&](double& x1, double& x2) {
    auto [i,j,k] = split_index_3d(v,f1.shape());
    x1 = 100*i + 10*j + k;
    x2 = 2*x1;
    v++;
  });

  auto time_function = [](const auto& function) -> double {
    runtime::trace::StopWatch stopwatch;
    for (size_t j=0; j<5; ++j) {
      function();
    }
    size_t N=10;
    stopwatch.start();
    for (size_t j=0; j<N; ++j) {
      function();
    }
    stopwatch.stop();
    return stopwatch.elapsed()/double(N);
  };

  double time_for_seq = time_function([&]{
    auto v1 = array::make_view<const double,3>(f1);
    auto v2 = array::make_view<const double,3>(f2);
    auto v3 = array::make_view<double,3>(f3);
    for (size_t i=0; i<ni; ++i) {
      for (size_t j=0; j<nj; ++j) {
        for (size_t k=0; k<nk; ++k) {
          v3(i,j,k) = v1(i,j,k) + v2(i,j,k);
        }
      }
    }
  });
  Log::info() << "timing with for loops seq  = " << time_for_seq << std::endl;

  double time_for_par = time_function([&]{
    auto v1 = array::make_view<const double,3>(f1);
    auto v2 = array::make_view<const double,3>(f2);
    auto v3 = array::make_view<double,3>(f3);
    atlas_omp_parallel_for (size_t i=0; i<ni; ++i) {
      for (size_t j=0; j<nj; ++j) {
        for (size_t k=0; k<nk; ++k) {
          v3(i,j,k) = v1(i,j,k) + v2(i,j,k);
        }
      }
    }
  });
  Log::info() << "timing with for loops par  = " << time_for_par << std::endl;

  double time_seq = time_function([&]{
    field::for_each_value(execution::seq, f1, f2, f3, [&](const double& x1, const double& x2, double& x3) {
      x3 = x1 + x2;
    });
  });
  Log::info() << "timing with execution::seq = " << time_seq << std::endl;

  double time_par = time_function([&]{
    field::for_each_value(execution::par, f1, f2, f3, [&](const double& x1, const double& x2, double& x3) {
      x3 = x1 + x2;
    });
  });
  Log::info() << "timing with execution::par = " << time_par << std::endl;
}



}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) { return atlas::test::run(argc, argv); }
