/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <map>
#include <string>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/output/Gmsh.h"
#include "atlas/output/detail/GmshIO.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

using atlas::Field;
using atlas::FieldSet;
using atlas::FunctionSpace;
using atlas::Mesh;
using eckit::Parametrisation;

namespace atlas {
namespace output {

// -----------------------------------------------------------------------------

std::string GmshFileStream::parallelPathName( const PathName& path, int part ) {
    std::stringstream s;
    // s << path.dirName() << "/" << path.baseName(false) << "_p" << part <<
    // ".msh";
    s << path.asString() << ".p" << part;
    return s.str();
}

// -----------------------------------------------------------------------------

GmshFileStream::GmshFileStream( const PathName& file_path, const char* mode, int part ) {
    PathName par_path( file_path );
    std::ios_base::openmode omode = std::ios_base::out;
    if ( std::string( mode ) == "w" ) {
        omode = std::ios_base::out;
    }
    else if ( std::string( mode ) == "a" ) {
        omode = std::ios_base::app;
    }

    if ( part < 0 || mpi::comm().size() == 1 ) {
        std::ofstream::open( file_path.localPath(), omode );
    }
    else {
        if ( mpi::comm().rank() == 0 ) {
            PathName par_path( file_path );
            std::ofstream par_file( par_path.localPath(), std::ios_base::out );
            for ( size_t p = 0; p < mpi::comm().size(); ++p ) {
                par_file << "Merge \"" << parallelPathName( file_path, p ) << "\";" << std::endl;
            }
            par_file.close();
        }
        PathName path( parallelPathName( file_path, part ) );
        std::ofstream::open( path.localPath(), omode );
    }
}

// -----------------------------------------------------------------------------

namespace detail {

// -----------------------------------------------------------------------------

void Gmsh::defaults() {
    config_.binary   = false;
    config_.nodes    = "xy";
    config_.gather   = false;
    config_.ghost    = false;
    config_.elements = true;
    config_.edges    = false;
    config_.levels.clear();
    config_.file        = "output.msh";
    config_.info        = false;
    config_.openmode    = "w";
    config_.coordinates = "xy";
}

// -----------------------------------------------------------------------------

namespace /*anonymous*/ {

// -----------------------------------------------------------------------------

void merge( Gmsh::Configuration& present, const eckit::Parametrisation& update ) {
    update.get( "binary", present.binary );
    update.get( "nodes", present.nodes );
    update.get( "gather", present.gather );
    update.get( "ghost", present.ghost );
    update.get( "elements", present.elements );
    update.get( "edges", present.edges );
    update.get( "levels", present.levels );
    update.get( "file", present.file );
    update.get( "info", present.info );
    update.get( "openmode", present.openmode );
    update.get( "coordinates", present.coordinates );
}

// -----------------------------------------------------------------------------

detail::GmshIO writer( const Gmsh::Configuration& c ) {
    detail::GmshIO gmsh;
    Gmsh::setGmshConfiguration( gmsh, c );
    return gmsh;
}

// -----------------------------------------------------------------------------

std::ios_base::openmode openmode( const Gmsh::Configuration& c ) {
    std::ios_base::openmode omode( std::ios_base::out );
    if ( std::string( c.openmode ) == "w" ) {
        omode = std::ios_base::out;
    }
    else if ( std::string( c.openmode ) == "a" ) {
        omode = std::ios_base::app;
    }
    if ( c.binary ) {
        omode |= std::ios::binary;
    }
    return omode;
}

// -----------------------------------------------------------------------------

}  // anonymous namespace

// -----------------------------------------------------------------------------

void Gmsh::setGmshConfiguration( detail::GmshIO& gmsh, const Gmsh::Configuration& c ) {
    gmsh.options.set( "ascii", not c.binary );
    gmsh.options.set( "nodes", c.nodes );
    gmsh.options.set( "gather", c.gather );
    gmsh.options.set( "ghost", c.ghost );
    gmsh.options.set( "elements", c.elements );
    gmsh.options.set( "edges", c.edges );
    gmsh.options.set( "levels", c.levels );
    gmsh.options.set( "info", c.info );
    gmsh.options.set( "nodes", c.coordinates );
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh( Stream& ) {
    defaults();
    ATLAS_NOTIMPLEMENTED;
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh( Stream&, const eckit::Parametrisation& config ) {
    defaults();
    merge( config_, config );
    ATLAS_NOTIMPLEMENTED;
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh( const PathName& file, const std::string& mode ) {
    defaults();
    config_.file     = file.asString();
    config_.openmode = std::string( mode );
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh( const PathName& file, const std::string& mode, const eckit::Parametrisation& config ) {
    defaults();
    merge( config_, config );
    config_.file     = file.asString();
    config_.openmode = std::string( mode );
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh( const PathName& file ) {
    defaults();
    config_.file = file.asString();
}

// -----------------------------------------------------------------------------

Gmsh::Gmsh( const PathName& file, const eckit::Parametrisation& config ) {
    defaults();
    merge( config_, config );
    config_.file = file.asString();
}

// -----------------------------------------------------------------------------

Gmsh::~Gmsh() {}

// -----------------------------------------------------------------------------

void Gmsh::write( const Mesh& mesh, const eckit::Parametrisation& config ) const {
    Gmsh::Configuration c = config_;
    merge( c, config );

    if ( c.coordinates == "xyz" and not mesh.nodes().has_field( "xyz" ) ) {
        Log::debug() << "Building xyz representation for nodes" << std::endl;
        mesh::actions::BuildXYZField( "xyz" )( const_cast<Mesh&>( mesh ) );
    }

    writer( c ).write( mesh, c.file );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void Gmsh::write( const Field& field, const eckit::Parametrisation& config ) const {
    Gmsh::Configuration c = config_;
    merge( c, config );
    writer( c ).write( field, c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void Gmsh::write( const FieldSet& fields, const eckit::Parametrisation& config ) const {
    Gmsh::Configuration c = config_;
    merge( c, config );
    writer( c ).write( fields, fields.field( 0 ).functionspace(), c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void Gmsh::write( const Field& field, const FunctionSpace& functionspace, const eckit::Parametrisation& config ) const {
    Gmsh::Configuration c = config_;
    merge( c, config );
    writer( c ).write( field, functionspace, c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void Gmsh::write( const FieldSet& fields, const FunctionSpace& functionspace,
                  const eckit::Parametrisation& config ) const {
    Gmsh::Configuration c = config_;
    merge( c, config );
    writer( c ).write( fields, functionspace, c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

extern "C" {

Gmsh* atlas__output__Gmsh__create_pathname_mode( const char* pathname, const char* mode ) {
    return new Gmsh( std::string( pathname ), std::string( mode ) );
}
Gmsh* atlas__output__Gmsh__create_pathname_mode_config( const char* pathname, const char* mode,
                                                        const Parametrisation* config ) {
    return new Gmsh( std::string( pathname ), std::string( mode ), *config );
}

}  // extern C

static OutputBuilder<detail::Gmsh> __gmsh( "gmsh" );

void force_link_atlas_output_detail_gmsh( void* ) {
    force_link_atlas_output_detail_gmsh( &__gmsh );
}

}  // namespace detail

//----------------------------------------------------------------------------------------------------------------------

Gmsh::Gmsh( const Output& output ) : Output( output ) {}

Gmsh::Gmsh( Stream& s ) : Output( new detail::Gmsh( s ) ) {}

Gmsh::Gmsh( Stream& s, const eckit::Parametrisation& c ) : Output( new detail::Gmsh( s, c ) ) {}

Gmsh::Gmsh( const PathName& p, const std::string& mode ) : Output( new detail::Gmsh( p, mode ) ) {}

Gmsh::Gmsh( const PathName& p, const std::string& mode, const eckit::Parametrisation& c ) :
    Output( new detail::Gmsh( p, mode, c ) ) {}

Gmsh::Gmsh( const PathName& p ) : Output( new detail::Gmsh( p ) ) {}

Gmsh::Gmsh( const PathName& p, const eckit::Parametrisation& c ) : Output( new detail::Gmsh( p, c ) ) {}

}  // namespace output
}  // namespace atlas
