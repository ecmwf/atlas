/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/detail/partitioner/Partitioner.h"

#include <map>
#include <string>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/util/Config.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/grid/detail/partitioner/CheckerboardPartitioner.h"
#include "atlas/grid/detail/partitioner/EqualRegionsPartitioner.h"
#include "atlas/grid/detail/partitioner/MatchingFunctionSpacePartitionerLonLatPolygon.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitioner.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerBruteForce.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerLonLatPolygon.h"
#include "atlas/grid/detail/partitioner/MatchingMeshPartitionerSphericalPolygon.h"
#include "atlas/library/config.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

#if ATLAS_HAVE_TRANS
#include "atlas/grid/detail/partitioner/TransPartitioner.h"
#endif

namespace {

static eckit::Mutex* local_mutex                                                       = nullptr;
static std::map<std::string, atlas::grid::detail::partitioner::PartitionerFactory*>* m = nullptr;
static pthread_once_t once                                                             = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, atlas::grid::detail::partitioner::PartitionerFactory*>();
}
}  // namespace

namespace atlas {
namespace grid {
namespace detail {
namespace partitioner {

Partitioner::Partitioner() : nb_partitions_( mpi::size() ) {}

Partitioner::Partitioner( const idx_t nb_partitions ) : nb_partitions_( nb_partitions ) {}

Partitioner::~Partitioner() = default;

idx_t Partitioner::nb_partitions() const {
    return nb_partitions_;
}

Distribution Partitioner::partition( const Grid& grid ) const {
    return Distribution( grid, atlas::grid::Partitioner( this ) );
}

namespace {

template <typename T>
void load_builder() {
    PartitionerBuilder<T>( "tmp" );
}

struct force_link {
    force_link() {
        load_builder<EqualRegionsPartitioner>();
        load_builder<CheckerboardPartitioner>();
#if ATLAS_HAVE_TRANS
        load_builder<TransPartitioner>();
#endif
    }
};

}  // namespace

PartitionerFactory::PartitionerFactory( const std::string& name ) : name_( name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    if ( m->find( name ) != m->end() ) {
        throw_Exception( "Partitioner with name [" + name + "] is already registered.", Here() );
    }
    ( *m )[name] = this;
}

PartitionerFactory::~PartitionerFactory() {
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    m->erase( name_ );
}

void PartitionerFactory::list( std::ostream& out ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    const char* sep = "";
    for ( std::map<std::string, PartitionerFactory*>::const_iterator j = m->begin(); j != m->end(); ++j ) {
        out << sep << ( *j ).first;
        sep = ", ";
    }
}

bool PartitionerFactory::has( const std::string& name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    return ( m->find( name ) != m->end() );
}

Partitioner* PartitionerFactory::build( const std::string& name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::map<std::string, PartitionerFactory*>::const_iterator j = m->find( name );

    Log::debug() << "Looking for PartitionerFactory [" << name << "]" << '\n';

    if ( j == m->end() ) {
        Log::error() << "No PartitionerFactory for [" << name << "]" << '\n';
        Log::error() << "PartitionerFactories are:" << '\n';
        for ( j = m->begin(); j != m->end(); ++j ) {
            Log::error() << "   " << ( *j ).first << '\n';
        }
        throw_Exception( std::string( "No PartitionerFactory called " ) + name );
    }

    return ( *j ).second->make();
}

Partitioner* PartitionerFactory::build( const std::string& name, const idx_t nb_partitions ) {
    atlas::util::Config config;
    return build (name, nb_partitions, config);
}

Partitioner* PartitionerFactory::build( const std::string& name, const idx_t nb_partitions, const eckit::Parametrisation & config ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::map<std::string, PartitionerFactory*>::const_iterator j = m->find( name );

    Log::debug() << "Looking for PartitionerFactory [" << name << "]" << '\n';

    if ( j == m->end() ) {
        Log::error() << "No PartitionerFactory for [" << name << "]" << '\n';
        Log::error() << "PartitionerFactories are:" << '\n';
        for ( j = m->begin(); j != m->end(); ++j ) {
            Log::error() << "   " << ( *j ).first << '\n';
        }
        throw_Exception( std::string( "No PartitionerFactory called " ) + name );
    }

    return ( *j ).second->make( nb_partitions, config );
}


Partitioner* MatchingPartitionerFactory::build( const std::string& type, const Mesh& partitioned ) {
    if ( type == MatchingMeshPartitionerSphericalPolygon::static_type() ) {
        return new MatchingMeshPartitionerSphericalPolygon( partitioned );
    }
    else if ( type == MatchingMeshPartitionerLonLatPolygon::static_type() ) {
        return new MatchingMeshPartitionerLonLatPolygon( partitioned );
    }
    else if ( type == MatchingMeshPartitionerBruteForce::static_type() ) {
        return new MatchingMeshPartitionerBruteForce( partitioned );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

Partitioner* MatchingPartitionerFactory::build( const std::string& type, const FunctionSpace& partitioned ) {
    if ( type == MatchingFunctionSpacePartitionerLonLatPolygon::static_type() ) {
        return new MatchingFunctionSpacePartitionerLonLatPolygon( partitioned );
    }
    else {
        ATLAS_NOTIMPLEMENTED;
    }
}

}  // namespace partitioner
}  // namespace detail
}  // namespace grid
}  // namespace atlas
