/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "TransInterface.h"

#include <cstring>

#include "atlas/array.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/Spectral.h"
#include "atlas/library/config.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/trans/Trans.h"
#include "atlas/trans/detail/TransFactory.h"
#include "atlas/trans/detail/TransImpl.h"

namespace atlas {
namespace trans {

class TransInterface {
private:
    const TransImpl* trans_;

public:
    TransInterface( const TransImpl* trans ) : trans_( trans ) {}
    int handle() const { return trans_->handle(); }
};

///////////////////////////////////////////////////////////////////////////////

extern "C" {

int atlas__Trans__has_backend( const char* backend ) {
    return Trans::hasBackend( std::string( backend ) );
}
void atlas__Trans__set_backend( const char* backend ) {
    Trans::backend( std::string( backend ) );
}
void atlas__Trans__backend( char*& backend, size_t& size ) {
    std::string s = Trans::backend();
    size          = s.size() + 1;
    backend       = new char[size];
    strcpy( backend, s.c_str() );
}


TransImpl* atlas__Trans__new( const Grid::Implementation* grid, int truncation ) {
    ATLAS_ASSERT( grid != nullptr, "Grid must not be null" );
    TransImpl* trans;
    {
        Grid g( grid );
        Trans t( g, truncation );
        trans = t.get();
        trans->attach();
    }
    trans->detach();
    return trans;
}

TransImpl* atlas__Trans__new_config( const Grid::Implementation* grid, int truncation,
                                     const eckit::Configuration* config ) {
    ATLAS_ASSERT( grid != nullptr, "Grid must not be null" );
    ATLAS_ASSERT( config != nullptr, "config must not be null" );
    TransImpl* trans;
    {
        Grid g( grid );
        Trans t( g, truncation, *config );
        trans = t.get();
        trans->attach();
    }
    trans->detach();
    return trans;
}


void atlas__Trans__delete( TransImpl* This ) {
    ATLAS_ASSERT( This != nullptr );
    delete This;
}

int atlas__Trans__handle( const TransImpl* This ) {
    ATLAS_ASSERT( This != nullptr );
    return TransInterface( This ).handle();
}


void atlas__Trans__invtrans_scalar( const TransImpl* t, int nb_fields, double scalar_spectra[],
                                    double scalar_fields[] ) {
    ATLAS_ASSERT( t != nullptr );
    return t->invtrans( nb_fields, scalar_spectra, scalar_fields );
}

void atlas__Trans__invtrans_vordiv2wind( const TransImpl* t, int nb_fields, double vorticity_spectra[],
                                         double divergence_spectra[], double wind_fields[] ) {
    ATLAS_ASSERT( t != nullptr );
    return t->invtrans( nb_fields, vorticity_spectra, divergence_spectra, wind_fields );
}

void atlas__Trans__dirtrans_scalar( const TransImpl* t, int nb_fields, double scalar_fields[],
                                    double scalar_spectra[] ) {
    ATLAS_ASSERT( t != nullptr );
    return t->dirtrans( nb_fields, scalar_fields, scalar_spectra );
}

void atlas__Trans__dirtrans_wind2vordiv( const TransImpl* t, int nb_fields, double wind_fields[],
                                         double vorticity_spectra[], double divergence_spectra[] ) {
    ATLAS_ASSERT( t != nullptr );
    return t->dirtrans( nb_fields, wind_fields, vorticity_spectra, divergence_spectra );
}

int atlas__Trans__truncation( const TransImpl* This ) {
    ATLAS_ASSERT( This != nullptr );
    return This->truncation();
}

const Grid::Implementation* atlas__Trans__grid( const TransImpl* This ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( This->grid() );
    return This->grid().get();
}

const functionspace::FunctionSpaceImpl* atlas__Trans__spectral( const TransImpl* This ) {
    ATLAS_ASSERT( This != nullptr );
    const auto spectral = This->spectral();
    ATLAS_ASSERT( spectral );
    return spectral.get();
}


void atlas__Trans__dirtrans_fieldset( const TransImpl* This, const field::FieldSetImpl* gpfields,
                                      field::FieldSetImpl* spfields, const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( gpfields );
    ATLAS_ASSERT( spfields );
    ATLAS_ASSERT( parameters );
    FieldSet fspfields( spfields );
    This->dirtrans( gpfields, fspfields, *parameters );
}

void atlas__Trans__dirtrans_field( const TransImpl* This, const field::FieldImpl* gpfield, field::FieldImpl* spfield,
                                   const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( spfield );
    ATLAS_ASSERT( gpfield );
    ATLAS_ASSERT( parameters );
    Field fspfield( spfield );
    This->dirtrans( gpfield, fspfield, *parameters );
}

void atlas__Trans__invtrans_fieldset( const TransImpl* This, const field::FieldSetImpl* spfields,
                                      field::FieldSetImpl* gpfields, const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( spfields );
    ATLAS_ASSERT( gpfields );
    ATLAS_ASSERT( parameters );
    FieldSet fgpfields( gpfields );
    This->invtrans( spfields, fgpfields, *parameters );
}

void atlas__Trans__invtrans_field( const TransImpl* This, const field::FieldImpl* spfield, field::FieldImpl* gpfield,
                                   const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( spfield );
    ATLAS_ASSERT( gpfield );
    ATLAS_ASSERT( parameters );
    Field fgpfield( gpfield );
    This->invtrans( spfield, fgpfield, *parameters );
}

void atlas__Trans__dirtrans_wind2vordiv_field( const TransImpl* This, const field::FieldImpl* gpwind,
                                               field::FieldImpl* spvor, field::FieldImpl* spdiv,
                                               const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( gpwind );
    ATLAS_ASSERT( spvor );
    ATLAS_ASSERT( spdiv );
    ATLAS_ASSERT( parameters );
    Field fspvor( spvor );
    Field fspdiv( spdiv );
    This->dirtrans_wind2vordiv( gpwind, fspvor, fspdiv, *parameters );
}

void atlas__Trans__invtrans_vordiv2wind_field( const TransImpl* This, const field::FieldImpl* spvor,
                                               const field::FieldImpl* spdiv, field::FieldImpl* gpwind,
                                               const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( spvor );
    ATLAS_ASSERT( spdiv );
    ATLAS_ASSERT( gpwind );
    ATLAS_ASSERT( parameters );
    Field fgpwind( gpwind );
    This->invtrans_vordiv2wind( spvor, spdiv, fgpwind, *parameters );
}

void atlas__Trans__invtrans( const TransImpl* This, int nb_scalar_fields, double scalar_spectra[], int nb_vordiv_fields,
                             double vorticity_spectra[], double divergence_spectra[], double gp_fields[],
                             const eckit::Configuration* parameters ) {
    ATLAS_ASSERT( This != nullptr );
    This->invtrans( nb_scalar_fields, scalar_spectra, nb_vordiv_fields, vorticity_spectra, divergence_spectra,
                    gp_fields, *parameters );
}

void atlas__Trans__invtrans_grad_field( const TransImpl* This, const field::FieldImpl* spfield,
                                        field::FieldImpl* gpfield, const eckit::Configuration* config ) {
    ATLAS_ASSERT( This != nullptr );
    ATLAS_ASSERT( spfield );
    ATLAS_ASSERT( gpfield );
    Field fgpfield( gpfield );
    This->invtrans_grad( spfield, fgpfield, *config );
}
}

}  // namespace trans
}  // namespace atlas
