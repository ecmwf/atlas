/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor
 * does it submit to any jurisdiction.
 */

#include <map>
#include <numeric>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"
#include "eckit/utils/Hash.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/HybridElements.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/DelaunayMeshGenerator.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/meshgenerator/StructuredMeshGenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/ErrorHandling.h"
#include "atlas/runtime/Log.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

namespace {

static eckit::Mutex* local_mutex                       = 0;
static std::map<std::string, MeshGeneratorFactory*>* m = 0;
static pthread_once_t once                             = PTHREAD_ONCE_INIT;

static void init() {
    local_mutex = new eckit::Mutex();
    m           = new std::map<std::string, MeshGeneratorFactory*>();
}

template <typename T>
void load_builder() {
    MeshGeneratorBuilder<T>( "tmp" );
}

struct force_link {
    force_link() {
        load_builder<meshgenerator::detail::StructuredMeshGenerator>();
        load_builder<meshgenerator::DelaunayMeshGenerator>();
    }
};

}  // namespace

//----------------------------------------------------------------------------------------------------------------------

MeshGeneratorImpl::MeshGeneratorImpl() {}

MeshGeneratorImpl::~MeshGeneratorImpl() {}

Mesh MeshGeneratorImpl::operator()( const Grid& grid ) const {
    Mesh mesh;
    generate( grid, mesh );
    return mesh;
}

Mesh MeshGeneratorImpl::operator()( const Grid& grid, const grid::Distribution& distribution ) const {
    Mesh mesh;
    generate( grid, distribution, mesh );
    return mesh;
}

Mesh MeshGeneratorImpl::generate( const Grid& grid ) const {
    Mesh mesh;
    generate( grid, mesh );
    return mesh;
}

Mesh MeshGeneratorImpl::generate( const Grid& grid, const grid::Distribution& distribution ) const {
    Mesh mesh;
    generate( grid, distribution, mesh );
    return mesh;
}

//----------------------------------------------------------------------------------------------------------------------

void MeshGeneratorImpl::generateGlobalElementNumbering( Mesh& mesh ) const {
    size_t loc_nb_elems = mesh.cells().size();
    std::vector<size_t> elem_counts( mpi::comm().size() );
    std::vector<int> elem_displs( mpi::comm().size() );

    ATLAS_TRACE_MPI( ALLGATHER ) { mpi::comm().allGather( loc_nb_elems, elem_counts.begin(), elem_counts.end() ); }

    elem_displs.at( 0 ) = 0;
    for ( size_t jpart = 1; jpart < mpi::comm().size(); ++jpart ) {
        elem_displs.at( jpart ) = elem_displs.at( jpart - 1 ) + elem_counts.at( jpart - 1 );
    }

    gidx_t gid = 1 + elem_displs.at( mpi::comm().rank() );

    array::ArrayView<gidx_t, 1> glb_idx = array::make_view<gidx_t, 1>( mesh.cells().global_index() );

    for ( size_t jelem = 0; jelem < mesh.cells().size(); ++jelem ) {
        glb_idx( jelem ) = gid++;
    }

    size_t max_glb_idx = std::accumulate( elem_counts.begin(), elem_counts.end(), size_t( 0 ) );

    mesh.cells().global_index().metadata().set( "human_readable", true );
    mesh.cells().global_index().metadata().set( "min", 1 );
    mesh.cells().global_index().metadata().set( "max", max_glb_idx );
}

void MeshGeneratorImpl::setProjection( Mesh& mesh, const Projection& p ) const {
    mesh.setProjection( p );
}

void MeshGeneratorImpl::setGrid( Mesh& mesh, const Grid& g, const grid::Distribution& d ) const {
    mesh.setGrid( g );
    mesh.metadata().set( "distribution", d.type() );
}

//----------------------------------------------------------------------------------------------------------------------

MeshGeneratorFactory::MeshGeneratorFactory( const std::string& name ) : name_( name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    ASSERT( m->find( name ) == m->end() );
    ( *m )[name] = this;
}

MeshGeneratorFactory::~MeshGeneratorFactory() {
    eckit::AutoLock<eckit::Mutex> lock( local_mutex );
    m->erase( name_ );
}

void MeshGeneratorFactory::list( std::ostream& out ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    const char* sep = "";
    for ( std::map<std::string, MeshGeneratorFactory*>::const_iterator j = m->begin(); j != m->end(); ++j ) {
        out << sep << ( *j ).first;
        sep = ", ";
    }
}

const MeshGenerator::Implementation* MeshGeneratorFactory::build( const std::string& name ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::map<std::string, MeshGeneratorFactory*>::const_iterator j = m->find( name );

    Log::debug() << "Looking for MeshGeneratorFactory [" << name << "]" << std::endl;

    if ( j == m->end() ) {
        Log::error() << "No MeshGeneratorFactory for [" << name << "]" << std::endl;
        Log::error() << "MeshGeneratorFactories are:" << std::endl;
        for ( j = m->begin(); j != m->end(); ++j )
            Log::error() << "   " << ( *j ).first << std::endl;
        throw eckit::SeriousBug( std::string( "No MeshGeneratorFactory called " ) + name );
    }

    return ( *j ).second->make();
}

const MeshGenerator::Implementation* MeshGeneratorFactory::build( const std::string& name,
                                                                  const eckit::Parametrisation& param ) {
    pthread_once( &once, init );

    eckit::AutoLock<eckit::Mutex> lock( local_mutex );

    static force_link static_linking;

    std::map<std::string, MeshGeneratorFactory*>::const_iterator j = m->find( name );

    Log::debug() << "Looking for MeshGeneratorFactory [" << name << "]" << std::endl;

    if ( j == m->end() ) {
        Log::error() << "No MeshGeneratorFactory for [" << name << "]" << std::endl;
        Log::error() << "MeshGeneratorFactories are:" << std::endl;
        for ( j = m->begin(); j != m->end(); ++j )
            Log::error() << "   " << ( *j ).first << std::endl;
        throw eckit::SeriousBug( std::string( "No MeshGeneratorFactory called " ) + name );
    }

    return ( *j ).second->make( param );
}

//----------------------------------------------------------------------------------------------------------------------

extern "C" {

void atlas__MeshGenerator__delete( MeshGenerator::Implementation* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); delete This; );
}

const MeshGenerator::Implementation* atlas__MeshGenerator__create_noconfig( const char* name ) {
    const MeshGenerator::Implementation* meshgenerator( 0 );
    ATLAS_ERROR_HANDLING( {
        MeshGenerator m( std::string{name} );
        meshgenerator = m.get();
        meshgenerator->attach();
    } meshgenerator->detach(); );
    return meshgenerator;
}

const MeshGenerator::Implementation* atlas__MeshGenerator__create( const char* name,
                                                                   const eckit::Parametrisation* params ) {
    const MeshGenerator::Implementation* meshgenerator( 0 );
    ATLAS_ERROR_HANDLING( ASSERT( params ); {
        MeshGenerator m( std::string( name ), *params );
        meshgenerator = m.get();
        meshgenerator->attach();
    } meshgenerator->detach(); );
    return meshgenerator;
}

Mesh::Implementation* atlas__MeshGenerator__generate__grid_griddist( const MeshGenerator::Implementation* This,
                                                                     const Grid::Implementation* grid,
                                                                     const grid::Distribution::impl_t* distribution ) {
    ATLAS_ERROR_HANDLING( Mesh::Implementation * m; {
        Mesh mesh = This->generate( Grid( grid ), grid::Distribution( distribution ) );
        mesh.get()->attach();
        m = mesh.get();
    } m->detach();
                          return m; );
    return nullptr;
}

Mesh::Implementation* atlas__MeshGenerator__generate__grid( const MeshGenerator::Implementation* This,
                                                            const Grid::Implementation* grid ) {
    ATLAS_ERROR_HANDLING( Mesh::Implementation * m; {
        Mesh mesh = This->generate( Grid( grid ) );
        ;
        mesh.get()->attach();
        m = mesh.get();
    } m->detach();
                          return m; );
    return nullptr;
}
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

MeshGenerator::MeshGenerator() : meshgenerator_( nullptr ) {}

MeshGenerator::MeshGenerator( const Implementation* meshgenerator ) : meshgenerator_( meshgenerator ) {}

MeshGenerator::MeshGenerator( const MeshGenerator& meshgenerator ) : meshgenerator_( meshgenerator.meshgenerator_ ) {}

MeshGenerator::MeshGenerator( const std::string& key, const eckit::Parametrisation& params ) :
    meshgenerator_( meshgenerator::MeshGeneratorFactory::build( key, params ) ) {}

void MeshGenerator::hash( eckit::Hash& h ) const {
    return meshgenerator_->hash( h );
}

Mesh MeshGenerator::generate( const Grid& g, const grid::Distribution& d ) const {
    return meshgenerator_->generate( g, d );
}

Mesh MeshGenerator::generate( const Grid& g ) const {
    return meshgenerator_->generate( g );
}

Mesh MeshGenerator::operator()( const Grid& g, const grid::Distribution& d ) const {
    return meshgenerator_->operator()( g, d );
}

Mesh MeshGenerator::operator()( const Grid& g ) const {
    return meshgenerator_->operator()( g );
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
