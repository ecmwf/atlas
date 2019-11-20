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

#include "eckit/eckit_version.h"
#if 10000 * ECKIT_MAJOR_VERSION + 100 * ECKIT_MINOR_VERSION < 10400
#include "eckit/parser/JSON.h"
#else
#include "eckit/log/JSON.h"
#endif

#include "eckit/config/Resource.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/log/Bytes.h"
#include "eckit/log/Log.h"
#include "eckit/types/FloatCompare.h"

#include "atlas/grid.h"
#include "atlas/grid/detail/grid/GridFactory.h"
#include "atlas/runtime/AtlasTool.h"

//----------------------------------------------------------------------------------------------------------------------

struct AtlasGrids : public atlas::AtlasTool {
    bool serial() override { return true; }
    int execute( const Args& args ) override;
    std::string briefDescription() override { return "Catalogue of available built-in grids"; }
    std::string usage() override { return name() + " GRID [OPTION]... [--help,-h]"; }
    std::string longDescription() override {
        return "Catalogue of available built-in grids\n"
               "\n"
               "       Browse catalogue of grids\n"
               "\n"
               "       GRID: unique identifier for grid \n"
               "           Example values: N80, F40, O24, L32\n";
    }

    AtlasGrids( int argc, char** argv ) : AtlasTool( argc, argv ) {
        add_option(
            new SimpleOption<bool>( "list", "List all grids. The names are possible values for the GRID argument" ) );
        add_option( new SimpleOption<bool>( "info", "List information about GRID" ) );
        add_option( new SimpleOption<bool>( "json", "Export json" ) );
        add_option( new SimpleOption<bool>( "rtable", "Export IFS rtable" ) );
        add_option( new SimpleOption<bool>( "check", "Check grid" ) );
        add_option( new SimpleOption<bool>( "check-uid", "Check grid uid required" ) );
        add_option( new SimpleOption<bool>( "check-boundingbox", "Check grid bounding_box(n,w,s,e) required" ) );
    }
};

//------------------------------------------------------------------------------------------------------

int AtlasGrids::execute( const Args& args ) {
    using namespace atlas;

    std::string key = args.count() ? args( 0 ) : "";

    bool info = false;
    args.get( "info", info );

    bool json = false;
    args.get( "json", json );

    bool rtable = false;
    args.get( "rtable", rtable );

    bool check = false;
    args.get( "check", check );

    bool check_uid = false;
    args.get( "check-uid", check_uid );

    bool check_bbox = false;
    args.get( "check-boundingbox", check_bbox );

    bool list = false;
    args.get( "list", list );

    bool do_run = list || ( !key.empty() && ( info || json || rtable || check ) );

    if ( !key.empty() && !do_run ) {
        Log::error() << "Option wrong or missing after '" << key << "'" << std::endl;
    }

    if ( list ) {
        Log::info() << "usage: atlas-grids GRID [OPTION]... [--help]\n" << std::endl;
        Log::info() << "Available grids:" << std::endl;
        for ( const auto& key : grid::GridFactory::keys() ) {
            Log::info() << "  -- " << key << std::endl;
        }
    }

    if ( !key.empty() ) {
        eckit::PathName path{key};
        Grid grid = path.exists() ? Grid( Grid::Spec{path} ) : Grid( key );

        if ( !grid ) {
            return failed();
        }

        if ( info ) {
            Log::info() << "Grid " << key << std::endl;
            Log::info() << "   name:                               " << grid.name() << std::endl;
            Log::info() << "   uid:                                " << grid.uid() << std::endl;
            if ( auto gaussian = GaussianGrid( grid ) ) {
                Log::info() << "   Gaussian N number:                  " << gaussian.N() << std::endl;
            }
            Log::info() << "   number of points:                   " << grid.size() << std::endl;


            size_t memsize = grid.size() * sizeof( double );

            Log::info() << "   memory footprint per field (dp):    " << eckit::Bytes( memsize ) << std::endl;

            if ( auto structuredgrid = StructuredGrid( grid ) ) {
                if ( not grid.projection() ) {
                    double deg, km;

                    Log::info() << "   number of latitudes (N-S):          " << structuredgrid.ny() << std::endl;
                    Log::info() << "   number of longitudes (max):         " << structuredgrid.nxmax() << std::endl;

                    deg = ( structuredgrid.y().front() - structuredgrid.y().back() ) / ( structuredgrid.ny() - 1 );
                    km  = deg * 40075. / 360.;
                    Log::info() << "   approximate resolution N-S:         " << std::setw( 10 ) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    deg = 360. / static_cast<double>( structuredgrid.nx( structuredgrid.ny() / 2 ) );
                    km  = deg * 40075. / 360.;
                    Log::info() << "   approximate resolution E-W equator: " << std::setw( 10 ) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    deg = 360. * std::cos( structuredgrid.y( structuredgrid.ny() / 4 ) * M_PI / 180. ) /
                          static_cast<double>( structuredgrid.nx( structuredgrid.ny() / 4 ) );
                    km = deg * 40075. / 360.;
                    Log::info() << "   approximate resolution E-W midlat:  " << std::setw( 10 ) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    deg = 360. * std::cos( structuredgrid.y().front() * M_PI / 180. ) /
                          static_cast<double>( structuredgrid.nx().front() );
                    km = deg * 40075. / 360.;


                    Log::info() << "   approximate resolution E-W pole:    " << std::setw( 10 ) << std::fixed << deg
                                << " deg   " << km << " km " << std::endl;

                    Log::info() << "   spectral truncation -- linear:      " << structuredgrid.ny() - 1 << std::endl;
                    Log::info() << "   spectral truncation -- quadratic:   "
                                << static_cast<int>( std::floor( 2. / 3. * structuredgrid.ny() + 0.5 ) ) - 1
                                << std::endl;
                    Log::info() << "   spectral truncation -- cubic:       "
                                << static_cast<int>( std::floor( 0.5 * structuredgrid.ny() + 0.5 ) ) - 1 << std::endl;
                }

                auto precision = Log::info().precision( 3 );
                if ( grid.projection().units() == "meters" ) {
                    Log::info() << "   x : [ " << std::setw( 10 ) << std::fixed << structuredgrid.xspace().min() / 1000.
                                << " , " << std::setw( 10 ) << std::fixed << structuredgrid.xspace().max() / 1000.
                                << " ] km" << std::endl;
                    Log::info() << "   y : [ " << std::setw( 10 ) << std::fixed << structuredgrid.yspace().min() / 1000.
                                << " , " << std::setw( 10 ) << std::fixed << structuredgrid.yspace().max() / 1000.
                                << " ] km" << std::endl;
                    if ( structuredgrid.xspace().nxmax() == structuredgrid.xspace().nxmin() ) {
                        Log::info() << "   dx : " << structuredgrid.xspace().dx()[0] / 1000. << " km" << std::endl;
                    }
                    Log::info() << "   dy : "
                                << std::abs( structuredgrid.yspace()[1] - structuredgrid.yspace()[0] ) / 1000. << " km"
                                << std::endl;
                    Log::info() << "   lonlat(centre)    : "
                                << grid.projection().lonlat(
                                       {0.5 * ( structuredgrid.xspace().max() + structuredgrid.xspace().min() ),
                                        0.5 * ( structuredgrid.yspace().max() + structuredgrid.yspace().min() )} )
                                << std::endl;
                    Log::info() << "   lonlat(xmin,ymax) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().min(), structuredgrid.yspace().max()} )
                                << std::endl;
                    Log::info() << "   lonlat(xmin,ymin) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().min(), structuredgrid.yspace().min()} )
                                << std::endl;
                    Log::info() << "   lonlat(xmax,ymin) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().max(), structuredgrid.yspace().min()} )
                                << std::endl;
                    Log::info() << "   lonlat(xmax,ymax) : "
                                << grid.projection().lonlat(
                                       {structuredgrid.xspace().max(), structuredgrid.yspace().max()} )
                                << std::endl;
                }
                if ( grid.projection().units() == "degrees" ) {
                    Log::info() << "   x : [ " << std::setw( 10 ) << std::fixed << structuredgrid.xspace().min()
                                << " , " << std::setw( 10 ) << std::fixed << structuredgrid.xspace().max() << " ] deg"
                                << std::endl;
                    Log::info() << "   y : [ " << std::setw( 10 ) << std::fixed << structuredgrid.yspace().min()
                                << " , " << std::setw( 10 ) << std::fixed << structuredgrid.yspace().max() << " ] deg"
                                << std::endl;
                }
                auto it = grid.lonlat().begin();
                Log::info() << "   lonlat(first)     : " << *it << std::endl;
                it += grid.size() - 1;
                Log::info() << "   lonlat(last)      : " << *it << std::endl;
                Log::info().precision( precision );
            }
        }

        if ( json ) {
            std::stringstream stream;
            eckit::JSON js( stream );
            js.precision( 16 );
            js << grid.spec();
            std::cout << stream.str() << std::endl;
        }

        if ( rtable ) {
            if ( auto structuredgrid = StructuredGrid( grid ) ) {
                std::stringstream stream;
                stream << "&NAMRGRI\n";
                for ( idx_t j = 0; j < structuredgrid.ny(); ++j ) {
                    stream << " NRGRI(" << std::setfill( '0' ) << std::setw( 5 ) << 1 + j << ")=" << std::setfill( ' ' )
                           << std::setw( 5 ) << structuredgrid.nx( j ) << ",\n";
                }
                stream << "/" << std::flush;
                std::cout << stream.str() << std::endl;
            }
        }

        if ( check ) {
            bool check_failed = false;
            Log::Channel out;
            out.setStream( Log::error() );

            eckit::PathName path{key};
            if ( not path.exists() ) {
                out << "Check failed:  " << key << " is not a file" << std::endl;
                return failed();
            }

            util::Config config_check;
            if ( not util::Config{path}.get( "check", config_check ) ) {
                out << "Check failed:  no \"check\" section in " << key << std::endl;
                return failed();
            }

            idx_t size;
            if ( config_check.get( "size", size ) ) {
                if ( grid.size() != size ) {
                    out << "Check failed: grid size " << grid.size() << " expected to be " << size << std::endl;
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for size skipped" << std::endl;
            }

            std::string uid;
            if ( config_check.get( "uid", uid ) ) {
                if ( grid.uid() != uid ) {
                    out << "Check failed: grid uid " << grid.uid() << " expected to be " << uid << std::endl;
                    check_failed = true;
                }
            }
            else if ( check_uid && uid.empty() ) {
                out << "Check failed: grid uid " << grid.uid() << " was not encoded in the check" << std::endl;
                check_failed = true;
            }
            else {
                Log::warning() << "Check for uid skipped" << std::endl;
            }


            auto equal = []( double a, double b ) { return eckit::types::is_approximately_equal( a, b, 5.e-4 ); };

            std::vector<double> first_point_lonlat;
            if ( config_check.get( "lonlat(first)", first_point_lonlat ) ) {
                PointLonLat first_point = *grid.lonlat().begin();
                if ( not equal( first_point.lon(), first_point_lonlat[0] ) or
                     not equal( first_point.lat(), first_point_lonlat[1] ) ) {
                    out << "Check failed: lonlat(first) " << first_point << " expected to be "
                        << PointLonLat( first_point_lonlat.data() ) << std::endl;
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for lonlat(first) skipped" << std::endl;
            }

            std::vector<double> last_point_lonlat;
            if ( config_check.get( "lonlat(last)", last_point_lonlat ) ) {
                PointLonLat last_point = *( grid.lonlat().begin() + ( grid.size() - 1 ) );
                if ( not equal( last_point.lon(), last_point_lonlat[0] ) or
                     not equal( last_point.lat(), last_point_lonlat[1] ) ) {
                    out << "Check failed: lonlat(last) " << last_point << " expected to be "
                        << PointLonLat( last_point_lonlat.data() ) << std::endl;
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for lonlat(last) skipped" << std::endl;
            }

            std::vector<double> bbox;
            if ( config_check.get( "bounding_box(n,w,s,e)", bbox ) && bbox.size() == 4 ) {
                auto bb = grid.lonlatBoundingBox();
                if ( !bb ) {
                    check_failed = true;
                    out << "Check failed: cannot calculate bounding box for " << grid.spec() << std::endl;
                }
                else if ( !equal( bb.north(), bbox[0] ) ) {
                    check_failed = true;
                    out << "Check failed: n=" << bb.north() << " expected to be " << bbox[0] << std::endl;
                }
                else if ( !equal( bb.west(), bbox[1] ) ) {
                    check_failed = true;
                    out << "Check failed: w=" << bb.west() << " expected to be " << bbox[1] << std::endl;
                }
                else if ( !equal( bb.south(), bbox[2] ) ) {
                    check_failed = true;
                    out << "Check failed: s=" << bb.south() << " expected to be " << bbox[2] << std::endl;
                }
                else if ( !equal( bb.east(), bbox[3] ) ) {
                    check_failed = true;
                    out << "Check failed: e=" << bb.east() << " expected to be " << bbox[3] << std::endl;
                }
            }
            else if ( check_bbox && bbox.size() != 4 ) {
                out << "Check failed: grid bounding_box(n,w,s,e) " << grid.lonlatBoundingBox()
                    << " was not encoded in the check" << std::endl;
                check_failed = true;
            }
            else {
                Log::warning() << "Check for bounding_box(n,w,s,e) skipped" << std::endl;
            }

            auto rel_equal = []( double a, double b ) { return std::abs( ( a - b ) / a ) < 1.e-6; };

            double xmin;
            if ( config_check.get( "xmin", xmin ) ) {
                if ( !rel_equal( RectangularDomain( grid.domain() ).xmin(), xmin ) ) {
                    auto precision = out.precision( 2 );
                    out << "Check failed: grid xmin " << std::fixed << RectangularDomain( grid.domain() ).xmin()
                        << " expected to be " << std::fixed << xmin << std::endl;
                    out.precision( precision );
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for xmin skipped" << std::endl;
            }

            double ymin;
            if ( config_check.get( "ymin", ymin ) ) {
                if ( !rel_equal( RectangularDomain( grid.domain() ).ymin(), ymin ) ) {
                    auto precision = out.precision( 2 );
                    out << "Check failed: grid ymin " << std::fixed << RectangularDomain( grid.domain() ).ymin()
                        << " expected to be " << std::fixed << ymin << std::endl;
                    out.precision( precision );
                    check_failed = true;
                }
            }
            else {
                Log::warning() << "Check for ymin skipped" << std::endl;
            }


            if ( check_failed ) {
                return failed();
            }
            Log::info() << "SUCCESS: All checks passed" << std::endl;
        }
    }
    return success();
}

//------------------------------------------------------------------------------------------------------

int main( int argc, char** argv ) {
    AtlasGrids tool( argc, argv );
    return tool.start();
}
