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
#include "atlas/output/detail/GmshIO.h"
#include "atlas/output/detail/GmshImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace output {
namespace detail {

// -----------------------------------------------------------------------------

void GmshImpl::defaults() {
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

    config_.configured_land_water = false;
    config_.land                  = false;
    config_.water                 = false;
}

// -----------------------------------------------------------------------------

namespace /*anonymous*/ {

// -----------------------------------------------------------------------------

void merge( GmshImpl::Configuration& present, const eckit::Parametrisation& update ) {
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
    if ( update.has( "water" ) || update.has( "land" ) ) {
        update.get( "land", present.land );
        update.get( "water", present.water );
        present.configured_land_water = true;
    }
}

// -----------------------------------------------------------------------------

detail::GmshIO writer( const GmshImpl::Configuration& c ) {
    detail::GmshIO gmsh;
    GmshImpl::setGmshConfiguration( gmsh, c );
    return gmsh;
}

// -----------------------------------------------------------------------------

std::ios_base::openmode openmode( const GmshImpl::Configuration& c ) {
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

void GmshImpl::setGmshConfiguration( detail::GmshIO& gmsh, const GmshImpl::Configuration& c ) {
    gmsh.options.set( "ascii", not c.binary );
    gmsh.options.set( "nodes", c.nodes );
    gmsh.options.set( "gather", c.gather );
    gmsh.options.set( "ghost", c.ghost );
    gmsh.options.set( "elements", c.elements );
    gmsh.options.set( "edges", c.edges );
    gmsh.options.set( "levels", c.levels );
    gmsh.options.set( "info", c.info );
    gmsh.options.set( "nodes", c.coordinates );
    if ( c.configured_land_water ) {
        gmsh.options.set( "land", c.land );
        gmsh.options.set( "water", c.water );
    }
}

// -----------------------------------------------------------------------------

GmshImpl::GmshImpl( std::ostream& ) {
    defaults();
    ATLAS_NOTIMPLEMENTED;  // JIRA issue ATLAS-254
}

// -----------------------------------------------------------------------------

GmshImpl::GmshImpl( std::ostream&, const eckit::Parametrisation& config ) {
    defaults();
    merge( config_, config );
    ATLAS_NOTIMPLEMENTED;  // JIRA issue ATLAS-254
}

// -----------------------------------------------------------------------------

GmshImpl::GmshImpl( const eckit::PathName& file, const std::string& mode ) {
    defaults();
    config_.file     = file.asString();
    config_.openmode = std::string( mode );
}

// -----------------------------------------------------------------------------

GmshImpl::GmshImpl( const eckit::PathName& file, const std::string& mode, const eckit::Parametrisation& config ) {
    defaults();
    merge( config_, config );
    config_.file     = file.asString();
    config_.openmode = std::string( mode );
}

// -----------------------------------------------------------------------------

GmshImpl::GmshImpl( const eckit::PathName& file ) {
    defaults();
    config_.file = file.asString();
}

// -----------------------------------------------------------------------------

GmshImpl::GmshImpl( const eckit::PathName& file, const eckit::Parametrisation& config ) {
    defaults();
    merge( config_, config );
    config_.file = file.asString();
}

// -----------------------------------------------------------------------------

GmshImpl::~GmshImpl() = default;

// -----------------------------------------------------------------------------

void GmshImpl::write( const Mesh& mesh, const eckit::Parametrisation& config ) const {
    GmshImpl::Configuration c = config_;
    merge( c, config );

    if ( c.coordinates == "xyz" and not mesh.nodes().has_field( "xyz" ) ) {
        Log::debug() << "Building xyz representation for nodes" << std::endl;
        mesh::actions::BuildXYZField( "xyz" )( const_cast<Mesh&>( mesh ) );
    }

    writer( c ).write( mesh, c.file );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void GmshImpl::write( const Field& field, const eckit::Parametrisation& config ) const {
    GmshImpl::Configuration c = config_;
    merge( c, config );
    writer( c ).write( field, c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void GmshImpl::write( const FieldSet& fields, const eckit::Parametrisation& config ) const {
    GmshImpl::Configuration c = config_;
    merge( c, config );
    writer( c ).write( fields, fields.field( 0 ).functionspace(), c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void GmshImpl::write( const Field& field, const FunctionSpace& functionspace,
                      const eckit::Parametrisation& config ) const {
    GmshImpl::Configuration c = config_;
    merge( c, config );
    writer( c ).write( field, functionspace, c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

void GmshImpl::write( const FieldSet& fields, const FunctionSpace& functionspace,
                      const eckit::Parametrisation& config ) const {
    GmshImpl::Configuration c = config_;
    merge( c, config );
    writer( c ).write( fields, functionspace, c.file, openmode( c ) );
    config_.openmode = "a";
}

// -----------------------------------------------------------------------------

static OutputBuilder<detail::GmshImpl> __gmsh( "gmsh" );

// -----------------------------------------------------------------------------

}  // namespace detail
}  // namespace output
}  // namespace atlas
