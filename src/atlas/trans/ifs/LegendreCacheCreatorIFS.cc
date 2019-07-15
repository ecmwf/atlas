/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/ifs/LegendreCacheCreatorIFS.h"
#include <sstream>
#include <string>
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/trans/Trans.h"
#include "eckit/utils/MD5.h"

namespace atlas {
namespace trans {

namespace {
static LegendreCacheCreatorBuilder<LegendreCacheCreatorIFS> builder( "ifs" );
}

namespace {

std::string truncate( const std::string& str ) {
    const int trunc = std::min( 10ul, str.size() );
    return str.substr( 0, trunc );
}

std::string hash( const Grid& grid ) {
    eckit::MD5 h;
    if ( StructuredGrid( grid ) && not grid.projection() ) {
        auto g = StructuredGrid( grid );
        h.add( g.y().data(), g.y().size() * sizeof( double ) );
    }
    else {
        grid.hash( h );
    }
    return truncate( h.digest() );
}

std::string hash( const eckit::Configuration& config ) {
    eckit::MD5 h;

    // Add options and other unique keys
    h << "flt" << config.getBool( "flt", false );

    return truncate( h.digest() );
}

}  // namespace

std::string LegendreCacheCreatorIFS::uid() const {
    if ( unique_identifier_.empty() ) {
        std::ostringstream stream;
        stream << "ifs-T" << truncation_ << "-";
        if ( GaussianGrid( grid_ ) ) {
            if ( RegularGaussianGrid( grid_ ) ) {
                stream << "RegularGaussianN" << GaussianGrid( grid_ ).N();
            }
            else {
                stream << "ReducedGaussianN" << GaussianGrid( grid_ ).N() << "-PL";
                stream << hash( grid_ );
            }
        }
        else if ( RegularLonLatGrid( grid_ ) ) {
            auto g = RegularLonLatGrid( grid_ );
            if ( g.standard() || g.shifted() ) {
                stream << ( g.standard() ? "L" : "S" ) << g.nx() << "x" << g.ny();
            }
            else {
                // We cannot make more assumptions on reusability for different grids
                stream << "grid-" << hash( grid_ );
            }
        }
        else {
            // We cannot make more assumptions on reusability for different grids
            stream << "grid-" << hash( grid_ );
        }
        stream << "-OPT" << hash( config_ );
        unique_identifier_ = stream.str();
    }
    return unique_identifier_;
}

LegendreCacheCreatorIFS::~LegendreCacheCreatorIFS() {}

bool LegendreCacheCreatorIFS::supported() const {
    if ( GaussianGrid( grid_ ) ) {
        return true;
    }
    else if ( RegularLonLatGrid( grid_ ) ) {
        auto g = RegularLonLatGrid( grid_ );
        if ( g.standard() || g.shifted() ) {
            return true;
        }
    }
    return false;
}

LegendreCacheCreatorIFS::LegendreCacheCreatorIFS( const Grid& grid, int truncation,
                                                  const eckit::Configuration& config ) :
    grid_( grid ),
    truncation_( truncation ),
    config_( config ) {}

void LegendreCacheCreatorIFS::create( const std::string& path ) const {
    Trans( grid_, truncation_, config_ | option::type( "ifs" ) | option::write_legendre( path ) );
}

Cache LegendreCacheCreatorIFS::create() const {
    return TransCache( Trans( grid_, truncation_, config_ | option::type( "ifs" ) ) );
}

size_t LegendreCacheCreatorIFS::estimate() const {
    return size_t( truncation_ * truncation_ * truncation_ ) / 2 * sizeof( double );
}


}  // namespace trans
}  // namespace atlas
