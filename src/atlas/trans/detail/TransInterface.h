/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/Grid.h"
#include "atlas/trans/detail/TransImpl.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace eckit {
class Configuration;
}  // namespace eckit

namespace atlas {
class Field;
class FieldSet;
namespace field {
class FieldImpl;
class FieldSetImpl;
}  // namespace field
namespace functionspace {
class FunctionSpaceImpl;
}
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

// C wrapper interfaces to C++ routines

extern "C" {

int atlas__Trans__has_backend(const char* backend);
void atlas__Trans__set_backend(const char* backend);
void atlas__Trans__backend(char*& backend, size_t& size);

TransImpl* atlas__Trans__new(const Grid::Implementation* grid, int truncation);
TransImpl* atlas__Trans__new_config(const Grid::Implementation* grid, int truncation,
                                    const eckit::Configuration* config);
void atlas__Trans__delete(TransImpl* trans);
void atlas__Trans__invtrans(const TransImpl* t, int nb_scalar_fields, double scalar_spectra[], int nb_vordiv_fields,
                            double vorticity_spectra[], double divergence_spectra[], double gp_fields[],
                            const eckit::Configuration* parameters);
void atlas__Trans__invtrans_scalar(const TransImpl* t, int nb_fields, double scalar_spectra[], double scalar_fields[]);
void atlas__Trans__invtrans_vordiv2wind(const TransImpl* t, int nb_fields, double vorticity_spectra[],
                                        double divergence_spectra[], double wind_fields[]);

void atlas__Trans__invtrans_adj(const TransImpl* t, int nb_scalar_fields, double gp_fields[], int nb_vordiv_fields,
                                double vorticity_spectra[], double divergence_spectra[], double scalar_spectra[],
                                const eckit::Configuration* parameters);
void atlas__Trans__invtrans_adj_scalar(const TransImpl* t, int nb_fields, double scalar_fields[],
                                       double scalar_spectra[]);
void atlas__Trans__invtrans_vordiv2wind_adj(const TransImpl* t, int nb_fields, double wind_fields[],
                                            double vorticity_spectra[], double divergence_spectra[]);

void atlas__Trans__dirtrans_scalar(const TransImpl* t, int nb_fields, double scalar_fields[], double scalar_spectra[]);
void atlas__Trans__dirtrans_wind2vordiv(const TransImpl* t, int nb_fields, double wind_fields[],
                                        double vorticity_spectra[], double divergence_spectra[]);
void atlas__Trans__dirtrans_wind2vordiv_field(const TransImpl* This, const field::FieldImpl* gpwind,
                                              field::FieldImpl* spvor, field::FieldImpl* spdiv,
                                              const eckit::Configuration* parameters);
// void atlas__Trans__specnorm( const TransImpl* t, int nb_fields, double spectra[], double norms[], int rank );
void atlas__Trans__dirtrans_fieldset(const TransImpl* This, const field::FieldSetImpl* gpfields,
                                     field::FieldSetImpl* spfields, const eckit::Configuration* parameters);
void atlas__Trans__dirtrans_field(const TransImpl* This, const field::FieldImpl* gpfield, field::FieldImpl* spfield,
                                  const eckit::Configuration* parameters);
void atlas__Trans__invtrans_fieldset(const TransImpl* This, const field::FieldSetImpl* spfields,
                                     field::FieldSetImpl* gpfields, const eckit::Configuration* parameters);
void atlas__Trans__invtrans_field(const TransImpl* This, const field::FieldImpl* spfield, field::FieldImpl* gpfield,
                                  const eckit::Configuration* parameters);
void atlas__Trans__invtrans_grad_field(const TransImpl* This, const field::FieldImpl* spfield,
                                       field::FieldImpl* gpfield, const eckit::Configuration* parameters);
void atlas__Trans__invtrans_vordiv2wind_field(const TransImpl* This, const field::FieldImpl* spvor,
                                              const field::FieldImpl* spdiv, field::FieldImpl* gpwind,
                                              const eckit::Configuration* parameters);
void atlas__Trans__invtrans_adj_fieldset(const TransImpl* This, const field::FieldSetImpl* gpfields,
                                         field::FieldSetImpl* spfields, const eckit::Configuration* parameters);
void atlas__Trans__invtrans_adj_field(const TransImpl* This, const field::FieldImpl* gpfield, field::FieldImpl* spfield,
                                      const eckit::Configuration* parameters);
void atlas__Trans__invtrans_grad_adj_field(const TransImpl* This, const field::FieldImpl* gpfield,
                                           field::FieldImpl* spfield, const eckit::Configuration* parameters);
void atlas__Trans__invtrans_vordiv2wind_adj_field(const TransImpl* This, const field::FieldImpl* gpwind,
                                                  field::FieldImpl* spvor, field::FieldImpl* spdiv,
                                                  const eckit::Configuration* parameters);

int atlas__Trans__handle(const TransImpl* trans);
int atlas__Trans__truncation(const TransImpl* This);
// int atlas__Trans__nspec2( const TransImpl* This );
// int atlas__Trans__nspec2g( const TransImpl* This );
// int atlas__Trans__ngptot( const TransImpl* This );
// int atlas__Trans__ngptotg( const TransImpl* This );
const Grid::Implementation* atlas__Trans__grid(const TransImpl* This);
const functionspace::FunctionSpaceImpl* atlas__Trans__spectral(const TransImpl* This);
}

// ------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
