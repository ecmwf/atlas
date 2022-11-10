/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


#include "eckit/config/Resource.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Bytes.h"
#include "eckit/log/JSON.h"
#include "eckit/log/Log.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/projection/detail/ProjectionFactory.h"
#include "atlas/runtime/AtlasTool.h"
#include "atlas/util/NormaliseLongitude.h"

namespace atlas {

template <typename Value>
class FixedFormat {
public:
    using value_type = Value;
    FixedFormat(value_type x, long precision): x_(x), precision_(precision > 0 ? precision : 20) {}
    void print(std::ostream& out) const {
        for (long precision = 0; precision <= precision_; ++precision) {
            if (is_precision(precision) || precision == precision_) {
                out << std::setprecision(precision);
                out << std::fixed << x_;
                break;
            }
        }
    }

    bool is_precision(long precision) const {
        std::stringstream ss;
        ss << std::setprecision(precision);
        ss << std::fixed << x_;
        value_type _x;
        ss >> _x;
        return std::abs(x_ - _x) < 1.e-20;
    }

    friend std::ostream& operator<<(std::ostream& out, const FixedFormat& This) {
        This.print(out);
        return out;
    }

private:
    float x_;
    long precision_;
};

FixedFormat<double> fixed_format(double x, long precision) {
    return FixedFormat<double>(x, precision);
}
FixedFormat<float> fixed_format(float x, long precision) {
    return FixedFormat<float>(x, precision);
}

}  // namespace atlas

//----------------------------------------------------------------------------------------------------------------------

struct AtlasGrids : public atlas::AtlasTool {
    bool serial() override { return true; }
    int execute(const Args& args) override;
    std::string briefDescription() override { return "Catalogue of available built-in grids"; }
    std::string usage() override { return name() + " <grid> [OPTION]... [--help,-h]"; }
    std::string longDescription() override {
        return "Catalogue of available built-in grids\n"
               "\n"
               "       Browse catalogue of grids\n"
               "\n"
               "       GRID: unique identifier for grid \n"
               "           Example values: N80, F40, O24, L32, CS-ED-12\n";
    }

    AtlasGrids(int argc, char** argv): AtlasTool(argc, argv) {
        add_option(
            new SimpleOption<bool>("list", "List all grids. The names are possible values for the <grid> argument"));
        add_option(new SimpleOption<bool>("info", "List information about <grid>"));
        add_option(new SimpleOption<bool>("json", "Export json"));
        add_option(new SimpleOption<bool>("rtable", "Export IFS rtable"));
        add_option(new SimpleOption<bool>("check", "Check grid"));
        add_option(new SimpleOption<bool>("check-uid", "Check grid uid required"));
        add_option(new SimpleOption<bool>("check-boundingbox", "Check grid bounding_box(n,w,s,e) required"));
        add_option(new SimpleOption<long>("precision", "Precision used for float output"));
        add_option(new SimpleOption<bool>("approximate-resolution", "Approximate resolution in degrees (North-South)"));
    }
};

//------------------------------------------------------------------------------------------------------

int AtlasGrids::execute(const Args& args) {
    using namespace atlas;

    std::string key = args.count() ? args(0) : "";

    bool info = false;
    args.get("info", info);

    bool json = false;
    args.get("json", json);

    bool rtable = false;
    args.get("rtable", rtable);

    bool check = false;
    args.get("check", check);

    bool check_uid = false;
    args.get("check-uid", check_uid);

    bool check_bbox = false;
    args.get("check-boundingbox", check_bbox);

    bool list = false;
    args.get("list", list);

    bool approximate_resolution = false;
    args.get("approximate-resolution", approximate_resolution);

    bool do_run = list || (!key.empty() && (info || json || rtable || check || approximate_resolution));

    if (!key.empty() && !do_run) {
        Log::error() << "Option wrong or missing after '" << key << "'" << std::endl;
    }

    if (list) {
        Log::info() << "usage: atlas-grids <grid> [OPTION]... [--help]\n" << std::endl;
        Log::info() << "\n";
        Log::info() << "Available grid types:" << std::endl;
        for (auto b : grid::GridBuilder::typeRegistry()) {
            Log::info() << "  -- " << b.second->type() << '\n';
        }
        Log::info() << "\n";
        Log::info() << "Available named grids:" << std::endl;

        size_t maxlen = 0;
        for (auto b : grid::GridBuilder::nameRegistry()) {
            for (const auto& name : b.second->names()) {
                maxlen = std::max(maxlen, name.size());
            }
        }

        for (auto b : grid::GridBuilder::nameRegistry()) {
            int c = 0;
            for (const auto& name : b.second->names()) {
                if (c == 0) {
                    Log::info() << "  -- " << std::left << std::setw(maxlen + 8) << name;
                    if (!b.second->type().empty()) {
                        Log::info() << "type: " << b.second->type();
                    }
                    Log::info() << std::endl;
                }
                else {
                    Log::info() << "     " << std::left << std::setw(maxlen + 8) << name << std::endl;
                }
                c++;
            }
        }
    }

    if (!key.empty()) {
        eckit::PathName path{key};
        Grid grid = path.exists() ? Grid(Grid::Spec{path}) : Grid(key);

        if (!grid) {
            return failed();
        }

        if (info) {
            Log::info() << "Grid " << key << std::endl;
            Log::info() << "   name:                               " << grid.name() << std::endl;
            Log::info() << "   uid:                                " << grid.uid() << std::endl;
            if (auto gaussian = GaussianGrid(grid)) {
                Log::info() << "   Gaussian N number:                  " << gaussian.N() << std::endl;
            }
            if (auto cubedsphere = CubedSphereGrid(grid)) {
                Log::info() << "   Cubedsphere faces:                   " << cubedsphere.N() << "x" << cubedsphere.N()
                            << "x6" << std::endl;
            }
            Log::info() << "   number of points:                   " << grid.size() << std::endl;

            Log::info() << "   memory footprint of grid:           " << eckit::Bytes(grid.footprint()) << std::endl;


            size_t memsize = grid.size() * sizeof(double);

            Log::info() << "   memory footprint per field (dp):    " << eckit::Bytes(memsize) << std::endl;

            if (auto structuredgrid = StructuredGrid(grid)) {
                if (not grid.projection()) {
                    double deg, km;

                    Log::info() << "   number of latitudes (N-S):          " << structuredgrid.ny() << std::endl;
                    Log::info() << "   number of longitudes (max):         " << structuredgrid.nxmax() << std::endl;

                    deg = (structuredgrid.y().front() - structuredgrid.y().back()) / (structuredgrid.ny() - 1);
                    km  = deg * 40075. / 360.;
                    Log::info() << "   approximate resolution N-S:         " << std::setw(10) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    deg = 360. / static_cast<double>(structuredgrid.nx(structuredgrid.ny() / 2));
                    km  = deg * 40075. / 360.;
                    Log::info() << "   approximate resolution E-W equator: " << std::setw(10) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    deg = 360. * std::cos(structuredgrid.y(structuredgrid.ny() / 4) * M_PI / 180.) /
                          static_cast<double>(structuredgrid.nx(structuredgrid.ny() / 4));
                    km = deg * 40075. / 360.;
                    Log::info() << "   approximate resolution E-W midlat:  " << std::setw(10) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    deg = 360. * std::cos(structuredgrid.y().front() * M_PI / 180.) /
                          static_cast<double>(structuredgrid.nx().front());
                    km = deg * 40075. / 360.;


                    Log::info() << "   approximate resolution E-W pole:    " << std::setw(10) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    Log::info() << "   spectral truncation -- linear:      " << structuredgrid.ny() - 1 << std::endl;
                    Log::info() << "   spectral truncation -- quadratic:   "
                                << static_cast<int>(std::floor(2. / 3. * structuredgrid.ny() + 0.5)) - 1 << std::endl;
                    Log::info() << "   spectral truncation -- cubic:       "
                                << static_cast<int>(std::floor(0.5 * structuredgrid.ny() + 0.5)) - 1 << std::endl;
                }

                auto precision = Log::info().precision(3);
                if (grid.projection().units() == "meters") {
                    Log::info() << "   x : [ " << std::setw(10) << std::fixed << structuredgrid.xspace().min() / 1000.
                                << " , " << std::setw(10) << std::fixed << structuredgrid.xspace().max() / 1000.
                                << " ] km" << std::endl;
                    Log::info() << "   y : [ " << std::setw(10) << std::fixed << structuredgrid.yspace().min() / 1000.
                                << " , " << std::setw(10) << std::fixed << structuredgrid.yspace().max() / 1000.
                                << " ] km" << std::endl;
                    if (structuredgrid.xspace().nxmax() == structuredgrid.xspace().nxmin()) {
                        Log::info() << "   dx : " << structuredgrid.xspace().dx()[0] / 1000. << " km" << std::endl;
                    }
                    Log::info() << "   dy : " << std::abs(structuredgrid.y(1) - structuredgrid.y(0)) / 1000. << " km"
                                << std::endl;
                    Log::info() << "   lonlat(centre)    : "
                                << grid.projection().lonlat(
                                       {0.5 * (structuredgrid.xspace().max() + structuredgrid.xspace().min()),
                                        0.5 * (structuredgrid.yspace().max() + structuredgrid.yspace().min())})
                                << std::endl;
                    Log::info() << "   lonlat(xmin,ymax) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().min(), structuredgrid.yspace().max()})
                                << std::endl;
                    Log::info() << "   lonlat(xmin,ymin) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().min(), structuredgrid.yspace().min()})
                                << std::endl;
                    Log::info() << "   lonlat(xmax,ymin) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().max(), structuredgrid.yspace().min()})
                                << std::endl;
                    Log::info() << "   lonlat(xmax,ymax) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().max(), structuredgrid.yspace().max()})
                                << std::endl;
                }
                if (grid.projection().units() == "degrees") {
                    Log::info() << "   x : [ " << std::setw(10) << std::fixed << structuredgrid.xspace().min() << " , "
                                << std::setw(10) << std::fixed << structuredgrid.xspace().max() << " ] deg"
                                << std::endl;
                    double ymin = std::min(structuredgrid.yspace().front(), structuredgrid.yspace().back());
                    double ymax = std::max(structuredgrid.yspace().front(), structuredgrid.yspace().back());
                    Log::info() << "   y : [ " << std::setw(10) << std::fixed << ymin << " , " << std::setw(10)
                                << std::fixed << ymax << " ] deg" << std::endl;
                }
                auto it = grid.lonlat().begin();
                Log::info() << "   lonlat(first)     : " << *it << std::endl;
                it += grid.size() - 1;
                Log::info() << "   lonlat(last)      : " << *it << std::endl;

                if (auto bb = grid.lonlatBoundingBox()) {
                    Log::info() << "   bounding_box(n,w,s,e) : { " << bb.north() << ", " << bb.west() << ", "
                                << bb.south() << ", " << bb.east() << " }" << std::endl;
                }

                Log::info().precision(precision);
            }
        }

        if (approximate_resolution) {
            if (auto structuredgrid = StructuredGrid(grid)) {
                if (structuredgrid.domain().global()) {
                    auto deg = (structuredgrid.y().front() - structuredgrid.y().back()) / (structuredgrid.ny() - 1);

                    long precision = -1;
                    args.get("precision", precision);
                    Log::info() << fixed_format(deg, precision) << std::endl;
                }
                else {
                    ATLAS_NOTIMPLEMENTED;
                }
            }
            else {
                ATLAS_NOTIMPLEMENTED;
            }
        }

        if (json) {
            std::stringstream stream;
            eckit::JSON js(stream);
            js.precision(16);
            js << grid.spec();
            std::cout << stream.str() << std::endl;
        }

        if (rtable) {
            if (auto structuredgrid = StructuredGrid(grid)) {
                std::stringstream stream;
                stream << "&NAMRGRI\n";
                for (idx_t j = 0; j < structuredgrid.ny(); ++j) {
                    stream << " NRGRI(" << std::setfill('0') << std::setw(5) << 1 + j << ")=" << std::setfill(' ')
                           << std::setw(5) << structuredgrid.nx(j) << ",\n";
                }
                stream << "/" << std::flush;
                std::cout << stream.str() << std::endl;
            }
        }

        if (check) {
            bool check_failed = false;
            Log::Channel out;
            out.setStream(Log::error());

            eckit::PathName path{key};
            if (not path.exists()) {
                out << "Check failed:  " << key << " is not a file" << std::endl;
                return failed();
            }

            util::Config config_check;
            if (not util::Config{path}.get("check", config_check)) {
                out << "Check failed:  no \"check\" section in " << key << std::endl;
                return failed();
            }

            bool strict = config_check.getBool("strict", true);

            idx_t size;
            if (config_check.get("size", size)) {
                if (grid.size() != size) {
                    out << "Check failed: grid size " << grid.size() << " expected to be " << size << std::endl;
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for size skipped" << std::endl;
            }

            std::string uid;
            if (config_check.get("uid", uid)) {
                if (uid == "ignore") {
                    out << "WARNING: ignoring uid explicitly" << std::endl;
                }
                else if (grid.uid() != uid) {
                    out << "Check failed: grid uid " << grid.uid() << " expected to be " << uid << std::endl;
                    check_failed = true;
                }
            }
            else if (check_uid && uid.empty()) {
                out << "Check failed: grid uid " << grid.uid() << " was not encoded in the check" << std::endl;
                check_failed = true;
            }
            else {
                Log::warning() << "Check for uid skipped" << std::endl;
            }


            auto equal       = [](double a, double b) { return eckit::types::is_approximately_equal(a, b, 5.e-4); };
            auto point_equal = [&equal](const PointLonLat& a, const PointLonLat& b) -> bool {
                return equal(a.lon(), b.lon()) && equal(a.lat(), b.lat());
            };
            auto difference_normalised = [](double a, double ref) {
                util::NormaliseLongitude normalised{ref - 180};
                return normalised(a) - ref;
            };
            auto point_equal_normalised = [&](const PointLonLat& a, const PointLonLat& b) -> bool {
                return equal(a.lat(), b.lat()) && equal(difference_normalised(a.lon(), b.lon()), 0.);
            };
            auto equal_normalised = [&](double a, double b) { return equal(difference_normalised(a, b), 0.); };

            auto check_lonlat = [&](const std::string& key, std::function<PointLonLat()>&& get_lonlat) {
                std::vector<double> lonlat_config;
                if (config_check.get(key, lonlat_config)) {
                    PointLonLat lonlat_check = {lonlat_config.data()};
                    PointLonLat lonlat       = get_lonlat();
                    if (not point_equal_normalised(lonlat, lonlat_check)) {
                        out << std::setprecision(4) << std::fixed << "Check failed: " << key << " " << lonlat
                            << " expected to be " << lonlat_check;
                        out << " ( normalised difference: {" << difference_normalised(lonlat.lon(), lonlat_check.lon())
                            << "," << lonlat.lat() - lonlat_check.lat() << "} )" << std::endl;
                        check_failed = true;
                    }
                    else if (strict && not point_equal(lonlat, lonlat_check)) {
                        out << std::setprecision(4) << std::fixed << "Check failed: " << key << " " << lonlat
                            << " expected to be " << lonlat_check;
                        out << " ( normalised difference: {" << difference_normalised(lonlat.lon(), lonlat_check.lon())
                            << "," << lonlat.lat() - lonlat_check.lat() << "} )" << std::endl;
                        check_failed = true;
                    }
                }
                else {
                    Log::warning() << "Check for " << key << " skipped" << std::endl;
                }
            };

            check_lonlat("lonlat(first)", [&]() { return grid.lonlat().front(); });
            check_lonlat("lonlat(last)", [&]() { return grid.lonlat().back(); });

            std::vector<double> bbox;
            if (config_check.get("bounding_box(n,w,s,e)", bbox) && bbox.size() == 4) {
                auto bb = grid.lonlatBoundingBox();
                if (!bb) {
                    check_failed = true;
                    out << "Check failed: cannot calculate bounding box for " << grid.spec() << std::endl;
                }
                else {
                    bool any_value_failed = false;
                    if (!equal(bb.north(), bbox[0])) {
                        any_value_failed = true;
                        out << "Check failed: n=" << bb.north() << " expected to be " << bbox[0] << std::endl;
                    }
                    if (!equal_normalised(bb.west(), bbox[1])) {
                        any_value_failed = true;
                        out << "Check failed: w=" << bb.west() << " expected to be " << bbox[1]
                            << " ( normalised difference : " << difference_normalised(bb.west(), bbox[1]) << " )"
                            << std::endl;
                    }
                    else if (strict && not equal(bb.west(), bbox[1])) {
                        out << "Check failed: w=" << bb.west() << " expected to be " << bbox[1]
                            << " ( normalised difference : " << difference_normalised(bb.west(), bbox[1]) << " )"
                            << std::endl;
                    }
                    if (!equal(bb.south(), bbox[2])) {
                        any_value_failed = true;
                        out << "Check failed: s=" << bb.south() << " expected to be " << bbox[2] << std::endl;
                    }
                    if (!equal_normalised(bb.east(), bbox[3])) {
                        any_value_failed = true;
                        out << "Check failed: e=" << bb.east() << " expected to be " << bbox[3]
                            << " ( normalised difference : " << difference_normalised(bb.east(), bbox[3]) << " )"
                            << std::endl;
                    }
                    else if (strict && not equal(bb.east(), bbox[3])) {
                        any_value_failed = true;
                        out << "Check failed: e=" << bb.east() << " expected to be " << bbox[3]
                            << " ( normalised difference : " << difference_normalised(bb.east(), bbox[3]) << " )"
                            << std::endl;
                    }
                    if (any_value_failed) {
                        check_failed = true;
                        out << "Check failed: bounding_box(n,w,s,e) [" << bb.north() << ", " << bb.west() << ", "
                            << bb.south() << ", " << bb.east() << "] expected to be [" << bbox[0] << ", " << bbox[1]
                            << ", " << bbox[2] << ", " << bbox[3] << "]" << std::endl;
                    }
                }
            }
            else if (check_bbox && bbox.size() != 4) {
                out << "Check failed: grid bounding_box(n,w,s,e) " << grid.lonlatBoundingBox()
                    << " was not encoded in the check" << std::endl;
                check_failed = true;
            }
            else {
                Log::warning() << "Check for bounding_box(n,w,s,e) skipped" << std::endl;
            }

            auto rel_equal = [](double a, double b) { return std::abs((a - b) / a) < 1.e-6; };

            double xmin;
            if (config_check.get("xmin", xmin)) {
                if (!rel_equal(RectangularDomain(grid.domain()).xmin(), xmin)) {
                    auto precision = out.precision(2);
                    out << "Check failed: grid xmin " << std::fixed << RectangularDomain(grid.domain()).xmin()
                        << " expected to be " << std::fixed << xmin << std::endl;
                    out.precision(precision);
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for xmin skipped" << std::endl;
            }

            double ymin;
            if (config_check.get("ymin", ymin)) {
                if (!rel_equal(RectangularDomain(grid.domain()).ymin(), ymin)) {
                    auto precision = out.precision(2);
                    out << "Check failed: grid ymin " << std::fixed << RectangularDomain(grid.domain()).ymin()
                        << " expected to be " << std::fixed << ymin << std::endl;
                    out.precision(precision);
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for ymin skipped" << std::endl;
            }

            if (projection::ProjectionFactory::has("proj")) {
                std::string proj_str;
                if (config_check.get("proj", proj_str)) {
                    Projection proj(util::Config("type", "proj") | util::Config("proj", proj_str));
                    auto check_proj = [&](const PointLonLat& lonlat, const PointLonLat& proj_lonlat,
                                          const std::string& point = "") {
                        if (not point_equal(lonlat, proj_lonlat)) {
                            if (point_equal_normalised(lonlat, proj_lonlat)) {
                                Log::warning()
                                    << "WARNING: Projection of " << point
                                    << " grid point is different from Proj only due to different normalisation: "
                                    << lonlat << " != " << proj_lonlat << std::endl;
                            }
                            else {
                                Log::info() << "Check failed: Projection of " << point
                                            << " grid point is different from Proj: " << lonlat << " != " << proj_lonlat
                                            << " normalised difference = {" << std::fixed << std::setprecision(12)
                                            << difference_normalised(lonlat.lon(), proj_lonlat.lon()) << ","
                                            << lonlat.lat() - proj_lonlat.lat() << "}" << std::endl;
                                check_failed = true;
                            }
                        }
                    };
                    check_proj(grid.lonlat().front(), proj.lonlat(grid.xy().front()), "first");
                    check_proj(grid.lonlat().back(), proj.lonlat(grid.xy().back()), "last");
                }
            }

            auto check_roundtrip = [&](const PointLonLat& point, bool strict = false) {
                using eckit::types::is_approximately_equal;
                const auto projection       = grid.projection();
                PointLonLat point_roundtrip = projection.lonlat(projection.xy(point));

                bool roundtrip = is_approximately_equal(point_roundtrip.lat(), point.lat(), 1.e-6) &&
                                 is_approximately_equal(point_roundtrip.lon(), point.lon(), 1.e-6);
                constexpr util::NormaliseLongitude normalised;
                bool normalised_roundtrip =
                    is_approximately_equal(point_roundtrip.lat(), point.lat(), 1.e-6) &&
                    is_approximately_equal(normalised(point_roundtrip.lon()), normalised(point.lon()));

                if (not normalised_roundtrip) {
                    check_failed = true;
                    Log::info() << "Check failed: Roundtrip of point " << point
                                << " failed, even with normalisation: " << point_roundtrip << " != " << point
                                << " normalised difference = {" << std::fixed << std::setprecision(12)
                                << difference_normalised(point_roundtrip.lon(), point.lon()) << ","
                                << point_roundtrip.lat() - point.lat() << "}" << std::endl;
                }
                if (strict && not roundtrip) {
                    check_failed = true;
                    Log::info() << "Check failed: Roundtrip of point " << point << " failed";
                    if (normalised_roundtrip) {
                        Log::info() << " but was OK with normalisation";
                    }
                    Log::info() << ": " << point_roundtrip << " != " << point << " normalised difference = {"
                                << std::fixed << std::setprecision(12)
                                << difference_normalised(point_roundtrip.lon(), point.lon()) << ","
                                << point_roundtrip.lat() - point.lat() << "}" << std::endl;
                }
            };

            std::vector<util::Config> roundtrip;
            if (config_check.get("roundtrip", roundtrip)) {
                for (auto& entry : roundtrip) {
                    std::vector<double> lonlat;
                    entry.get("lonlat", lonlat);
                    check_roundtrip({lonlat[0], lonlat[1]}, strict);
                }
            }


            if (check_failed) {
                return failed();
            }
            Log::info() << "SUCCESS: All checks passed" << std::endl;
        }
    }
    return success();
}

//------------------------------------------------------------------------------------------------------

int main(int argc, char** argv) {
    AtlasGrids tool(argc, argv);
    return tool.start();
}
