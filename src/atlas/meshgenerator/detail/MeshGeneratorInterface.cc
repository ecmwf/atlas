/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/meshgenerator/detail/MeshGeneratorInterface.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"
#include "atlas/runtime/ErrorHandling.h"

using atlas::Mesh;

namespace atlas {
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

extern "C" {

void atlas__MeshGenerator__delete( MeshGenerator::Implementation* This ) {
    ATLAS_ERROR_HANDLING( ASSERT( This ); delete This; );
}

const MeshGenerator::Implementation* atlas__MeshGenerator__create_noconfig( const char* name ) {
    const MeshGenerator::Implementation* meshgenerator( nullptr );
    ATLAS_ERROR_HANDLING( {
        MeshGenerator m( std::string{name} );
        meshgenerator = m.get();
        meshgenerator->attach();
    } meshgenerator->detach(); );
    return meshgenerator;
}

const MeshGenerator::Implementation* atlas__MeshGenerator__create( const char* name,
                                                                   const eckit::Parametrisation* params ) {
    const MeshGenerator::Implementation* meshgenerator( nullptr );
    ATLAS_ERROR_HANDLING( ASSERT( params ); {
        MeshGenerator m( std::string( name ), *params );
        meshgenerator = m.get();
        meshgenerator->attach();
    } meshgenerator->detach(); );
    return meshgenerator;
}

Mesh::Implementation* atlas__MeshGenerator__generate__grid_griddist(
    const MeshGenerator::Implementation* This, const Grid::Implementation* grid,
    const grid::Distribution::Implementation* distribution ) {
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
}  // extern "C"

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
