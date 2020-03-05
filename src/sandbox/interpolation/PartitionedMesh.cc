/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "PartitionedMesh.h"

#include "atlas/grid/Partitioner.h"
#include "atlas/output/Gmsh.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"

namespace atlas {
namespace interpolation {

PartitionedMesh::PartitionedMesh( const std::string& partitioner, const std::string& generator,
                                  bool generatorTriangulate, double generatorAngle ) :
    optionPartitioner_( partitioner ), optionGenerator_( generator ) {
    generatorParams_.set( "three_dimensional", false );
    generatorParams_.set( "patch_pole", true );
    generatorParams_.set( "include_pole", false );
    generatorParams_.set( "triangulate", generatorTriangulate );
    generatorParams_.set( "angle", generatorAngle );
}

void PartitionedMesh::writeGmsh( const std::string& fileName, const FieldSet& fields ) {
    util::Config output_config;
    //output_config.set( "coordinates", std::string( "xyz" ) );
    output_config.set( "ghost", true );

    output::Gmsh out( fileName, output_config );
    out.write( mesh_ );

    if ( not fields.empty() ) {
        out.write( fields );
    }
}

void PartitionedMesh::partition( const Grid& grid ) {
    ATLAS_TRACE( "PartitionedMesh::partition()" );

    partitioner_ = Partitioner( optionPartitioner_ );

    MeshGenerator meshgen( optionGenerator_, generatorParams_ );
    mesh_ = meshgen.generate( grid, partitioner_.partition( grid ) );
}

void PartitionedMesh::partition( const Grid& grid, const PartitionedMesh& other ) {
    ATLAS_TRACE( "PartitionedMesh::partition(other)" );

    partitioner_ = grid::MatchingMeshPartitioner( other.mesh_, util::Config( "type", optionPartitioner_ ) );

    MeshGenerator meshgen( optionGenerator_, generatorParams_ );
    mesh_ = meshgen.generate( grid, partitioner_.partition( grid ) );
}

}  // namespace interpolation
}  // namespace atlas
