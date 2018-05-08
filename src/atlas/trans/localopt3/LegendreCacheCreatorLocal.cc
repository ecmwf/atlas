/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/localopt3/LegendreCacheCreatorLocal.h"
#include <string>
#include <sstream>
#include "eckit/utils/MD5.h"
#include "atlas/grid.h"
#include "atlas/option.h"
#include "atlas/trans/Trans.h"

namespace atlas {
namespace trans {

namespace {
static LegendreCacheCreatorBuilder<LegendreCacheCreatorLocal> builder( "local" );
}

namespace {

std::string truncate( const std::string& str ) {
  const int trunc = std::min(10ul,str.size());
  return str.substr( 0, trunc );
}

std::string hash( const Grid& grid ) {
  eckit::MD5 h;
  if( grid::StructuredGrid( grid ) && not grid.projection() ) {
    auto g = grid::StructuredGrid( grid );
    h.add( g.y().data(), g.y().size() * sizeof(double) );
  } else {
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

}

std::string LegendreCacheCreatorLocal::uid() const {
  if( unique_identifier_.empty() ) {
    std::ostringstream stream;
    stream << "local-T" << truncation_ << "-";
    if( grid::GaussianGrid( grid_ ) ) {
      // Same cache for any global Gaussian grid
      stream << "GaussianN" << grid::GaussianGrid( grid_ ).N();
    } else if( grid::RegularLonLatGrid( grid_ ) ) {
      // Same cache for any global regular grid
      auto g = grid::RegularLonLatGrid( grid_ );
      stream << ( g.shiftedLat() ? "S" : "L" ) << "+x" << g.ny();
      // The above '+' is a placeholder for any g.nx()
    } else {
      // We cannot make more assumptions on reusability for different grids
      stream << "grid-" << hash( grid_ );
    }
    stream << "-OPT" << hash( config_ );
    unique_identifier_ = stream.str();
  }
  return unique_identifier_;
}

LegendreCacheCreatorLocal::~LegendreCacheCreatorLocal() {}

LegendreCacheCreatorLocal::LegendreCacheCreatorLocal( const Grid& grid, int truncation, const eckit::Configuration& config ) :
  grid_(grid),
  truncation_(truncation),
  config_(config) {
}

bool LegendreCacheCreatorLocal::supported() const {
  return true;
}

void LegendreCacheCreatorLocal::create( const std::string& path ) const {
  Trans( grid_, truncation_, config_ | option::type("local") | option::write_legendre( path ) );
}

Cache LegendreCacheCreatorLocal::create() const {
  NOTIMP;
}

size_t LegendreCacheCreatorLocal::estimate() const {
    return size_t(truncation_ * truncation_ * truncation_) / 2 * sizeof(double);
}


}  // namespace trans
}  // namespace atlas
