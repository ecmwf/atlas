/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "Grid.h"

#include <vector>

#include "eckit/memory/Factory.h"
#include "eckit/utils/MD5.h"

#include "atlas/grid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

static void checkSizeOfPoint() {
    // compile time check support C++11
    static_assert( sizeof( PointXY ) == 2 * sizeof( double ), "Grid requires size of Point to be 2*double" );

    // runtime check
    ASSERT( sizeof( PointXY ) == 2 * sizeof( double ) );
}

std::string Grid::className() {
    return "atlas.Grid";
}

const Grid* Grid::create( const Config& config ) {
    std::string name;
    if ( config.get( "name", name ) ) { return create( name, config ); }

    std::string type;
    if ( config.get( "type", type ) ) {
        const GridBuilder::Registry& registry = GridBuilder::typeRegistry();
        if ( registry.find( type ) != registry.end() ) {
            const GridBuilder& gc = *registry.at( type );
            return gc.create( config );
        }
    }

    if ( name.size() ) { Log::info() << "name provided: " << name << std::endl; }
    if ( type.size() ) { Log::info() << "type provided: " << type << std::endl; }
    if ( name.empty() && type.empty() ) { throw eckit::BadParameter( "no name or type in configuration", Here() ); }
    else {
        throw eckit::BadParameter( "name or type in configuration don't exist", Here() );
    }
}

const Grid* Grid::create( const std::string& name, const Grid::Config& config ) {
    const GridBuilder::Registry& registry = GridBuilder::nameRegistry();
    for ( GridBuilder::Registry::const_iterator it = registry.begin(); it != registry.end(); ++it ) {
        const Grid* grid = it->second->create( name, config );
        if ( grid ) { return grid; }
    }

    // Throw exception
    std::ostringstream log;
    log << "Could not construct Grid from the name \"" << name << "\"\n";
    log << "Accepted names are: \n";
    for ( GridBuilder::Registry::const_iterator it = registry.begin(); it != registry.end(); ++it ) {
        log << "  -  " << *it->second << "\n";
    }
    throw eckit::BadParameter( log.str() );
    //    return GridBuilder::createNamed(name);
}

const Grid* Grid::create( const Grid& grid, const Domain& domain ) {
    if ( grid.type() == "structured" ) {
        const Structured& g = dynamic_cast<const Structured&>( grid );
        return new Structured( g.name(), g.xspace(), g.yspace(), g.projection(), domain );
    }
    else {
        NOTIMP;
    }
}


Grid::Grid() {
    checkSizeOfPoint();
}

Grid::~Grid() {}

Grid::uid_t Grid::uid() const {
    if ( uid_.empty() ) { uid_ = hash(); }
    return uid_;
}

std::string Grid::hash() const {
    if ( hash_.empty() ) {
        eckit::MD5 md5;
        hash( md5 );
        hash_ = md5.digest();
    }
    return hash_;
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
