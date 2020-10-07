/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <string>

#include "atlas/grid/Distribution.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/runtime/Exception.h"

#include "atlas/meshgenerator/detail/MeshGeneratorFactory.h"
#include "atlas/meshgenerator/detail/MeshGeneratorImpl.h"

using atlas::Mesh;

namespace atlas {

//----------------------------------------------------------------------------------------------------------------------

MeshGenerator::MeshGenerator( const std::string& key, const eckit::Parametrisation& config ) :
    Handle( meshgenerator::MeshGeneratorFactory::build( key, config ) ) {}

MeshGenerator::MeshGenerator( const eckit::Parametrisation& config ) :
    Handle( meshgenerator::MeshGeneratorFactory::build(
        [&config]() {
            std::string key;
            ATLAS_ASSERT( config.get( "type", key ), "type must be specified in MeshGenerator configuration" );
            return key;
        }(),
        config ) ) {}

void MeshGenerator::hash( eckit::Hash& h ) const {
    return get()->hash( h );
}

Mesh MeshGenerator::generate( const Grid& g, const grid::Distribution& d ) const {
    return get()->generate( g, d );
}
Mesh MeshGenerator::generate( const Grid& g, const grid::Partitioner& p ) const {
    return get()->generate( g, p );
}
Mesh MeshGenerator::generate( const Grid& g ) const {
    return get()->generate( g );
}

Mesh MeshGenerator::operator()( const Grid& g, const grid::Distribution& d ) const {
    return get()->operator()( g, d );
}
Mesh MeshGenerator::operator()( const Grid& g, const grid::Partitioner& p ) const {
    return get()->operator()( g, p );
}
Mesh MeshGenerator::operator()( const Grid& g ) const {
    return get()->operator()( g );
}

std::string MeshGenerator::type() const {
    return get()->type();
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace atlas
