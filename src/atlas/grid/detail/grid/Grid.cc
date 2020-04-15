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

#include "eckit/utils/MD5.h"

#include "atlas/grid.h"
#include "atlas/grid/detail/grid/GridBuilder.h"
#include "atlas/grid/detail/grid/Structured.h"
#include "atlas/grid/detail/grid/Unstructured.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace grid {
namespace detail {
namespace grid {

static void checkSizeOfPoint() {
    // compile time check support C++11
    static_assert( sizeof( PointXY ) == 2 * sizeof( double ), "Grid requires size of Point to be 2*double" );

    // runtime check
    ATLAS_ASSERT( sizeof( PointXY ) == 2 * sizeof( double ) );
}

const Grid* Grid::create( const Config& config ) {
    std::string name;
    if ( config.get( "name", name ) ) {
        return create( name, config );
    }

    std::string type;
    if ( config.get( "type", type ) ) {
        const GridBuilder::Registry& registry = GridBuilder::typeRegistry();
        if ( registry.find( type ) != registry.end() ) {
            const GridBuilder& gc = *registry.at( type );
            return gc.create( config );
        }
    }

    if ( name.size() ) {
        Log::info() << "name provided: " << name << std::endl;
    }
    if ( type.size() ) {
        Log::info() << "type provided: " << type << std::endl;
    }
    if ( name.empty() && type.empty() ) {
        throw_Exception( "no name or type in configuration", Here() );
    }
    else {
        throw_Exception( "name or type in configuration don't exist", Here() );
    }
}

const Grid* Grid::create( const std::string& name ) {
    return create( name, util::NoConfig() );
}

const Grid* Grid::create( const std::string& name, const Grid::Config& config ) {
    const GridBuilder::Registry& registry = GridBuilder::nameRegistry();
    for ( GridBuilder::Registry::const_iterator it = registry.begin(); it != registry.end(); ++it ) {
        const Grid* grid = it->second->create( name, config );
        if ( grid ) {
            return grid;
        }
    }

    // Throw exception
    std::ostringstream log;
    log << "Could not construct Grid from the name \"" << name << "\"\n";
    log << "Accepted names are: \n";
    for ( GridBuilder::Registry::const_iterator it = registry.begin(); it != registry.end(); ++it ) {
        log << "  -  " << *it->second << "\n";
    }
    throw_Exception( log.str() );
    //    return GridBuilder::createNamed(name);
}

const Grid* Grid::create( const Grid& grid, const Domain& domain ) {
    if ( grid.type() == "structured" ) {
        const Structured& g = dynamic_cast<const Structured&>( grid );
        return new Structured( g.name(), g.xspace(), g.yspace(), g.projection(), domain );
    }
    else {
        return new Unstructured( grid, domain );
    }
}


Grid::Grid() {
    checkSizeOfPoint();
}

Grid::~Grid() {
    while ( grid_observers_.size() ) {
        GridObserver* o = grid_observers_.back();
        o->onGridDestruction( *this );
        o->unregisterGrid( *this );  // will also delete observer from mesh
    }
}

Grid::uid_t Grid::uid() const {
    if ( uid_.empty() ) {
        uid_ = hash();
    }
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

void Grid::attachObserver( GridObserver& observer ) const {
    if ( std::find( grid_observers_.begin(), grid_observers_.end(), &observer ) == grid_observers_.end() ) {
        grid_observers_.push_back( &observer );
    }
}

void Grid::detachObserver( GridObserver& observer ) const {
    grid_observers_.erase( std::remove( grid_observers_.begin(), grid_observers_.end(), &observer ),
                           grid_observers_.end() );
}

}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
