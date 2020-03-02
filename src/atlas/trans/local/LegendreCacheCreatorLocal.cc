/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/local/LegendreCacheCreatorLocal.h"

#include <cmath>
#include <sstream>
#include <string>

#include "eckit/types/FloatCompare.h"
#include "eckit/utils/MD5.h"

#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/runtime/Exception.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/local/TransLocal.h"

namespace atlas {
namespace trans {

namespace {
static LegendreCacheCreatorBuilder<LegendreCacheCreatorLocal> builder( "local" );
}

namespace {

std::string truncate( const std::string& str ) {
    const int trunc = std::min( 10ul, str.size() );
    return str.substr( 0, trunc );
}

std::string hash( const Grid& grid ) {
    eckit::MD5 h;

    StructuredGrid structured( grid );
    if ( structured && not grid.projection() ) {
        for ( auto& y : structured.y() ) {
            h.add( std::lround( y * 1.e8 ) );
        }
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

std::string LegendreCacheCreatorLocal::uid() const {
    if ( unique_identifier_.empty() ) {
        std::ostringstream stream;
        auto give_up = [&]() {
            // We cannot make more assumptions on reusability for different grids
            stream << "grid-" << hash( grid_ );
        };
        stream << "local-T" << truncation_ << "-";
        StructuredGrid structured( grid_ );
        if ( grid_.projection() ) {
            give_up();
        }
        else if ( GaussianGrid( grid_ ) ) {
            // Same cache for any global Gaussian grid
            stream << "GaussianN" << GaussianGrid( grid_ ).N();
        }
        else if ( RegularLonLatGrid( grid_ ) ) {
            // Same cache for any global regular grid
            auto g = RegularLonLatGrid( grid_ );

            const double dy_2 = 90. / double( g.ny() );
            bool shifted_lat  = eckit::types::is_approximately_equal( g.y().front(), 90. - dy_2 ) &&
                               eckit::types::is_approximately_equal( g.y().back(), -90. + dy_2 );
            bool standard_lat = eckit::types::is_approximately_equal( g.y().front(), 90. ) &&
                                eckit::types::is_approximately_equal( g.y().back(), -90. );

            if ( standard_lat ) {
                stream << "L"
                       << "-ny" << g.ny();
            }
            else if ( shifted_lat ) {
                stream << "S"
                       << "-ny" << g.ny();
            }
            else {  // I don't think we get here, but just in case, give up
                give_up();
            }
        }
        else if ( RegularGrid( grid_ ) && structured.yspace().type() == "linear" ) {
            RectangularDomain domain( grid_.domain() );
            ATLAS_ASSERT( domain );
            stream << "Regional";
            stream << "-south" << domain.ymin();
            stream << "-north" << domain.ymax();
            stream << "-ny" << structured.ny();
        }
        else {  // It gets too complicated, so let's not be smart
            give_up();
        }
        stream << "-OPT" << hash( config_ );
        unique_identifier_ = stream.str();
    }
    return unique_identifier_;
}

LegendreCacheCreatorLocal::~LegendreCacheCreatorLocal() = default;

LegendreCacheCreatorLocal::LegendreCacheCreatorLocal( const Grid& grid, int truncation,
                                                      const eckit::Configuration& config ) :
    grid_( grid ),
    truncation_( truncation ),
    config_( config ) {}

bool LegendreCacheCreatorLocal::supported() const {
    if ( not StructuredGrid( grid_ ) ) {
        return false;
    }
    if ( grid_.projection() ) {
        return false;
    }
    return true;
}

void LegendreCacheCreatorLocal::create( const std::string& path ) const {
    Trans tmp( grid_, truncation_, config_ | option::type( "local" ) | option::write_legendre( path ) );
}

Cache LegendreCacheCreatorLocal::create() const {
    util::Config export_legendre( "export_legendre", true );
    Trans tmp( grid_, truncation_, config_ | option::type( "local" ) | export_legendre );
    auto impl = dynamic_cast<const TransLocal*>( tmp.get() );
    return impl->export_legendre_;
}

size_t LegendreCacheCreatorLocal::estimate() const {
    return size_t( truncation_ * truncation_ * truncation_ ) / 2 * sizeof( double );
}


}  // namespace trans
}  // namespace atlas
