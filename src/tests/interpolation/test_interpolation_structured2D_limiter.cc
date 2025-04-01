/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <iomanip>

#include "atlas/library.h"
#include "atlas/runtime/Log.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/output/Gmsh.h"

#include "tests/AtlasTestEnvironment.h"

using atlas::functionspace::NodeColumns;
using atlas::functionspace::StructuredColumns;
using atlas::functionspace::PointCloud;
using atlas::util::Config;

namespace atlas {
namespace test {

//-----------------------------------------------------------------------------


std::vector<double> read_array(const std::string& filename) {
    std::string files_dir="/Users/willem/work/debug-ifs-mgrids-remap/files/H15p3_cubicatlasinterp/";
    std::ifstream file (files_dir+filename);
    assert(file.is_open());
    std::size_t size;
    file >> size;
    std::vector<double> array(size);
    for(std::size_t j=0; j<size; ++j) {
        file >> array[j];
    }
    return array;
}
std::vector<atlas::PointXY> to_points(const std::vector<double>& lon, const std::vector<double>& lat) {
    std::vector<atlas::PointXY> points(lon.size());
    for( std::size_t j=0; j<points.size(); ++j) {
        points[j][0] = lon[j];
        points[j][1] = lat[j];
    }
    return points;
}


std::set<atlas::idx_t> ignore_overshoot{6271, 6272, 7196, 7212};
std::set<atlas::idx_t> ignore_undershoot{43138,43319};

template <typename View>
std::pair<double,atlas::idx_t> view_max_loc(View v) {
    double m = -std::numeric_limits<double>::max();
    atlas::idx_t loc = 0;
    for(atlas::idx_t j=0; j<v.size(); ++j) {
        if(v[j] > m) {
            if (ignore_overshoot.find(j) == ignore_overshoot.end()) {
               m = v[j];
               loc = j;
            }
        }
    }
    return {m,loc};
}

std::pair<double,atlas::idx_t> view_max_loc(atlas::array::ArrayView<const double,2> v) {
    auto all = atlas::array::Range::all();
    return view_max_loc( v.slice(all,0) );
}

std::pair<double,atlas::idx_t> view_max_loc(atlas::array::ArrayView<const double,3> v) {
    auto all = atlas::array::Range::all();
    return view_max_loc( v.slice(all,0,0) );
}

std::pair<double,atlas::idx_t> max_loc(const Field& field) {
    if( field.rank() == 1 ) {
        return view_max_loc(atlas::array::make_view<double,1>(field));
    }
    else if( field.rank() == 2 ) {
        return view_max_loc(atlas::array::make_view<double,2>(field));
    }
    else if( field.rank() == 3 ) {
        return view_max_loc(atlas::array::make_view<double,3>(field));
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

template <typename View>
std::pair<double,atlas::idx_t> view_min_loc(View v) {
    double m = std::numeric_limits<double>::max();
    atlas::idx_t loc = 0;
    for(atlas::idx_t j=0; j<v.size(); ++j) {
        if(v[j] < m) {
            if (ignore_undershoot.find(j) == ignore_undershoot.end()) {
                m = v[j];
                loc = j;
            }
        }
    }
    return {m,loc};
}

std::pair<double,atlas::idx_t> view_min_loc(atlas::array::ArrayView<const double,2> v) {
    auto all = atlas::array::Range::all();
    return view_min_loc( v.slice(all,0) );
}

std::pair<double,atlas::idx_t> view_min_loc(atlas::array::ArrayView<const double,3> v) {
    auto all = atlas::array::Range::all();
    return view_min_loc( v.slice(all,0,0) );
}


std::pair<double,atlas::idx_t> min_loc(const Field& field) {
    if( field.rank() == 1 ) {
        return view_min_loc(atlas::array::make_view<double,1>(field));
    }
    else if( field.rank() == 2 ) {
        return view_min_loc(atlas::array::make_view<double,2>(field));
    }
    else if( field.rank() == 3 ) {
        return view_min_loc(atlas::array::make_view<double,3>(field));
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}


std::ostream& operator<<(std::ostream& out, const std::pair<double,atlas::idx_t>& mloc) {
    out << std::setw(13) << std::fixed << std::setprecision(10) << mloc.first << "  at index " << mloc.second;
    return out;
}

std::vector<atlas::PointXY> overshooting_dg_points() {
    return std::vector<atlas::PointXY>{
        {180., 45.2048},
        {0.,-88.9603}, // undershoot
        {360.,-88.9603} // undershoot
    };
}

std::vector<atlas::PointXY> read_dg_points() {
    static auto dg_lon = read_array("lambda_dg_nod.txt");
    static auto dg_lat = read_array("theta_dg_nod.txt");
    return to_points(dg_lon,dg_lat);
};

void read_tracer(atlas::Field& field_ifs) {
    // the tracer read from file should match closely the analytical function below 'init_tracer'
    static auto ifs_trc = read_array("trc_ifs_gp.txt");
    if( field_ifs.rank() == 1 ) {
        auto view_ifs = array::make_view<double,1>(field_ifs);
        for( std::size_t j=0; j<std::min<std::size_t>(ifs_trc.size(), field_ifs.shape(0)); ++j) {
            view_ifs(j) = ifs_trc[j];
        }
    }
    else if( field_ifs.rank() == 2 ) {
        auto view_ifs = array::make_view<double,2>(field_ifs);
        for( std::size_t j=0; j<std::min<std::size_t>(ifs_trc.size(), field_ifs.shape(0)); ++j) {
            view_ifs(j,0) = ifs_trc[j];
        }
    }
    else if( field_ifs.rank() == 3 ) {
        auto view_ifs = array::make_view<double,3>(field_ifs);
        for( std::size_t j=0; j<std::min<std::size_t>(ifs_trc.size(), field_ifs.shape(0)); ++j) {
            view_ifs(j,0,0) = ifs_trc[j];
        }
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    field_ifs.set_dirty();
}

void init_tracer(atlas::Field& field_ifs) {
    auto func = [](double lon, double lat) -> double {
        constexpr double lon_ctr = M_PI;
        constexpr double lat_ctr = M_PI/4.;
        constexpr double sigma_lon = 2.*M_PI/5.;
        constexpr double sigma_lat = 2.*M_PI/5.;
        constexpr double to_rad = M_PI/180.;
        lon *= to_rad;
        lat *= to_rad;
        auto sqr = [](double x){ return x*x;};
        return std::exp( - sqr( (lon-lon_ctr) / (std::sqrt(2.)*sigma_lon)) )
            * std::exp( - sqr( (lat-lat_ctr) / (std::sqrt(2.)*sigma_lat)) );
    };


    auto ifs_lonlat = atlas::array::make_view<double,2>(field_ifs.functionspace().lonlat());
    if( field_ifs.rank() == 1 ) {
        auto view_ifs = array::make_view<double,1>(field_ifs);
        for( std::size_t j=0; j<view_ifs.shape(0); ++j) {
            view_ifs(j) = func(ifs_lonlat(j,0), ifs_lonlat(j,1));
        }
    }
    else if( field_ifs.rank() == 2 ) {
        auto view_ifs = array::make_view<double,2>(field_ifs);
        for( std::size_t j=0; j<view_ifs.shape(0); ++j) {
            view_ifs(j,0) = func(ifs_lonlat(j,0), ifs_lonlat(j,1));
        }
    }
    else if( field_ifs.rank() == 3 ) {
        auto view_ifs = array::make_view<double,3>(field_ifs);
        for( std::size_t j=0; j<view_ifs.shape(0); ++j) {
            view_ifs(j,0,0) = func(ifs_lonlat(j,0), ifs_lonlat(j,1));
        }
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    field_ifs.set_dirty();
}

void init_crazy(atlas::Field& field) {
    if( field.rank() == 1 ) {
        array::make_view<double,1>(field).assign(99999999999.);
    }
    else if( field.rank() == 2 ) {
        array::make_view<double,2>(field).assign(99999999999.);
    }
    else if( field.rank() == 3 ) {
        array::make_view<double,3>(field).assign(99999999999.);
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}



// ----------------------------------------------------------------------

void do_test(bool limiter, int rank) {

    // std::vector<atlas::PointXY> dg_points = read_dg_points();
    //     This reads coordinates from files in ascii with format "nb_points value[0] value[1] ... value[nb_points-1]"
    std::vector<atlas::PointXY> dg_points = overshooting_dg_points();

    atlas::functionspace::PointCloud        dg_fs(dg_points);
    atlas::functionspace::StructuredColumns ifs_fs(atlas::Grid("O80"), atlas::option::halo(3));
    atlas::Field field_ifs;
    atlas::Field field_dg;

    if( rank == 1 ) {
        field_ifs = ifs_fs.createField<double>();
        field_dg  = dg_fs.createField<double>();
    }
    else if( rank == 2 ) {
        int nlev = 60;
        field_ifs = ifs_fs.createField<double>(atlas::option::levels(nlev));
        field_dg  = dg_fs.createField<double>(atlas::option::levels(nlev));
    }
    else if( rank == 3 ) {
        int nlev = 60;
        int ntrc = 1;
        field_ifs = ifs_fs.createField<double>(atlas::option::levels(nlev)|atlas::option::variables(ntrc));
        field_dg  = dg_fs.createField<double>(atlas::option::levels(nlev)|atlas::option::variables(ntrc));
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
    init_crazy(field_dg);
    init_tracer(field_ifs);
    // read_tracer(field_ifs);

    atlas::util::Config config;
    config.set("type", "structured-quasibicubic");
    config.set("limiter", limiter);
    config.set("matrix_free", true);
    atlas::Interpolation interpolator(config, ifs_fs, dg_fs);
    interpolator.execute(field_ifs, field_dg);

    auto max_ifs = max_loc(field_ifs);
    auto min_ifs = min_loc(field_ifs);
    auto max_dg  = max_loc(field_dg);
    auto min_dg  = min_loc(field_dg);

    int precision = 10;
    std::cout << "max_ifs:  " << max_ifs << std::endl;
    std::cout << "max_dg:   " << max_dg  << std::endl;
    std::cout << "max_diff: " << std::setw(13) << std::fixed << std::setprecision(precision) << max_ifs.first - max_dg.first << std::endl;
    std::cout << std::endl;
    std::cout << "min_ifs:  " << min_ifs << std::endl;
    std::cout << "min_dg:   " << min_dg  << std::endl;
    std::cout << "min_diff: " << std::setw(13) << std::fixed << std::setprecision(precision) << min_dg.first - min_ifs.first << std::endl;

    if( limiter && max_dg.first > max_ifs.first ) {
        std::cout << std::endl;
        std::cout << "FAILURE: overshoot in interpolation detected " << std::endl;
        auto lonlat = atlas::array::make_view<double,2>(dg_fs.lonlat());
        std::cout << " coordinate ["<<max_dg.second<<"]: " << atlas::PointLonLat{lonlat(max_dg.second,0),lonlat(max_dg.second,1)} << std::endl;
        EXPECT(max_dg.first <= max_ifs.first);
    }
    if( limiter && min_dg.first < min_ifs.first ) {
        std::cout << std::endl;
        std::cout << "FAILURE: undershoot in interpolation detected " << std::endl;
        auto lonlat = atlas::array::make_view<double,2>(dg_fs.lonlat());
        std::cout << " coordinate ["<<min_dg.second<<"]: " << atlas::PointLonLat{lonlat(min_dg.second,0),lonlat(min_dg.second,1)} << std::endl;
        EXPECT(min_dg.first >= min_ifs.first);
    }
}

CASE("structured-quasibicubic without limiter, rank 1") {
    bool limiter = false;
    int rank = 1;
    do_test(limiter, rank);
}
CASE("structured-quasibicubic with limiter, rank 1") {
    bool limiter = true;
    int rank = 1;
    do_test(limiter, rank);
}
CASE("structured-quasibicubic without limiter, rank 2") {
    bool limiter = false;
    int rank = 2;
    do_test(limiter, rank);
}
CASE("structured-quasibicubic with limiter, rank 2") {
    bool limiter = true;
    int rank = 2;
    do_test(limiter, rank);
}
CASE("structured-quasibicubic without limiter, rank 3") {
    bool limiter = false;
    int rank = 3;
    do_test(limiter, rank);
}
CASE("structured-quasibicubic with limiter, rank 3") {
    bool limiter = true;
    int rank = 3;
    do_test(limiter, rank);
}

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
