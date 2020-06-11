/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Partitioner.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/grid/detail/distribution/DistributionImpl.h"
#include "atlas/grid/detail/partitioner/Partitioner.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/option.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Trace.h"
#include "atlas/util/Config.h"

namespace atlas {
namespace grid {

using Factory = detail::partitioner::PartitionerFactory;

bool Partitioner::exists( const std::string& type ) {
    return Factory::has( type );
}

Partitioner::Partitioner( const std::string& type ) : Handle( Factory::build( type ) ) {}

Partitioner::Partitioner( const std::string& type, const idx_t nb_partitions ) :
    Handle( Factory::build( type, nb_partitions ) ) {}

namespace {
detail::partitioner::Partitioner* partitioner_from_config( const std::string& type,
                                                           const Partitioner::Config& config ) {
    long partitions = mpi::size();
    config.get( "partitions", partitions );
    return Factory::build( type, partitions, config );
}
detail::partitioner::Partitioner* partitioner_from_config( const Partitioner::Config& config ) {
    std::string type;
    if ( not config.get( "type", type ) ) {
        throw_Exception( "'type' missing in configuration for Partitioner", Here() );
    }
    return partitioner_from_config( type, config );
}
}  // namespace

Partitioner::Partitioner( const std::string& type, const Config& config ) :
    Handle( partitioner_from_config( type, config ) ) {}

Partitioner::Partitioner( const Config& config ) : Handle( partitioner_from_config( config ) ) {}

void Partitioner::partition( const Grid& grid, int part[] ) const {
    ATLAS_TRACE( "Partitioner::partition" );
    get()->partition( grid, part );
}

Distribution Partitioner::partition( const Grid& grid ) const {
    return get()->partition( grid );
}

idx_t Partitioner::nb_partitions() const {
    return get()->nb_partitions();
}

std::string Partitioner::type() const {
    return get()->type();
}

MatchingPartitioner::MatchingPartitioner() : Partitioner() {}

grid::detail::partitioner::Partitioner* matching_mesh_partititioner( const Mesh& mesh,
                                                                     const Partitioner::Config& config ) {
    std::string type( "lonlat-polygon" );
    config.get( "type", type );
    return grid::detail::partitioner::MatchingPartitionerFactory::build( type, mesh );
}

MatchingPartitioner::MatchingPartitioner( const Mesh& mesh ) : MatchingPartitioner( mesh, util::NoConfig() ) {}


MatchingPartitioner::MatchingPartitioner( const Mesh& mesh, const Config& config ) :
    Partitioner( matching_mesh_partititioner( mesh, config ) ) {}


grid::detail::partitioner::Partitioner* matching_functionspace_partititioner( const FunctionSpace& functionspace,
                                                                              const Partitioner::Config& config ) {
    std::string type( "lonlat-polygon" );
    config.get( "type", type );
    return grid::detail::partitioner::MatchingPartitionerFactory::build( type, functionspace );
}

MatchingPartitioner::MatchingPartitioner( const FunctionSpace& functionspace ) :
    MatchingPartitioner( functionspace, util::NoConfig() ) {}


MatchingPartitioner::MatchingPartitioner( const FunctionSpace& functionspace, const Config& config ) :
    Partitioner( matching_functionspace_partititioner( functionspace, config ) ) {}


extern "C" {

detail::partitioner::Partitioner* atlas__grid__Partitioner__new( const Partitioner::Config* config ) {
    detail::partitioner::Partitioner* p;
    {
        Partitioner partitioner( *config );
        p = const_cast<detail::partitioner::Partitioner*>( partitioner.get() );
        p->attach();
    }
    p->detach();
    return p;
}

detail::partitioner::Partitioner* atlas__grid__Partitioner__new_type( const char* type ) {
    detail::partitioner::Partitioner* p;
    {
        Partitioner partitioner{option::type( type )};
        p = const_cast<detail::partitioner::Partitioner*>( partitioner.get() );
        p->attach();
    }
    p->detach();
    return p;
}


detail::partitioner::Partitioner* atlas__grid__MatchingMeshPartitioner__new( const Mesh::Implementation* mesh,
                                                                             const Partitioner::Config* config ) {
    detail::partitioner::Partitioner* p;
    {
        MatchingPartitioner partitioner( Mesh( mesh ), *config );
        p = const_cast<detail::partitioner::Partitioner*>( partitioner.get() );
        p->attach();
    }
    p->detach();
    return p;
}

detail::partitioner::Partitioner* atlas__grid__MatchingFunctionSpacePartitioner__new(
    const FunctionSpace::Implementation* functionspace, const Partitioner::Config* config ) {
    detail::partitioner::Partitioner* p;
    {
        MatchingPartitioner partitioner( FunctionSpace( functionspace ), *config );
        p = const_cast<detail::partitioner::Partitioner*>( partitioner.get() );
        p->attach();
    }
    p->detach();
    return p;
}

Distribution::Implementation* atlas__grid__Partitioner__partition( const Partitioner::Implementation* This,
                                                                   const Grid::Implementation* grid ) {
    Distribution::Implementation* d;
    {
        Distribution distribution = This->partition( Grid( grid ) );
        d                         = const_cast<Distribution::Implementation*>( distribution.get() );
        d->attach();
    }
    d->detach();
    return d;
}

void atlas__grid__Partitioner__delete( detail::partitioner::Partitioner* This ) {
    delete This;
}

}  // extern "C"

}  // namespace grid
}  // namespace atlas
