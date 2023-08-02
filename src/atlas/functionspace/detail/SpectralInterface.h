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

#include "atlas/functionspace/Spectral.h"

// -------------------------------------------------------------------
// C wrapper interfaces to C++ routines
namespace atlas {
namespace field {
class FieldSetImpl;
class FieldImpl;
}  // namespace field
namespace trans {
class TransImpl;
}
namespace functionspace {

extern "C" {
const detail::Spectral* atlas__SpectralFunctionSpace__new__config(const eckit::Configuration* config);
void atlas__SpectralFunctionSpace__delete(detail::Spectral* This);
field::FieldImpl* atlas__fs__Spectral__create_field(const detail::Spectral* This, const eckit::Configuration* options);
void atlas__SpectralFunctionSpace__gather(const detail::Spectral* This, const field::FieldImpl* local,
                                          field::FieldImpl* global);
void atlas__SpectralFunctionSpace__gather_fieldset(const detail::Spectral* This, const field::FieldSetImpl* local,
                                                   field::FieldSetImpl* global);
void atlas__SpectralFunctionSpace__scatter(const detail::Spectral* This, const field::FieldImpl* global,
                                           field::FieldImpl* local);
void atlas__SpectralFunctionSpace__scatter_fieldset(const detail::Spectral* This, const field::FieldSetImpl* global,
                                                    field::FieldSetImpl* local);
void atlas__SpectralFunctionSpace__norm(const detail::Spectral* This, const field::FieldImpl* field, double norm[],
                                        int rank);
void atlas__SpectralFunctionSpace__nspec2(const detail::Spectral* This, int& nspec2);
void atlas__SpectralFunctionSpace__nspec2g(const detail::Spectral* This, int& nspec2g);
void atlas__SpectralFunctionSpace__truncation(const detail::Spectral* This, int& truncation);
void atlas__SpectralFunctionSpace__nump(const detail::Spectral* This, int& nump);
void atlas__SpectralFunctionSpace__nmyms(const detail::Spectral* This, const int*& nmyms, int& size);
void atlas__SpectralFunctionSpace__nasm0(const detail::Spectral* This, const int*& nasm0, int& size);
void atlas__SpectralFunctionSpace__nvalue(const detail::Spectral* This, const int*& nvalue, int& size);
void atlas__SpectralFunctionSpace__levels(const detail::Spectral* This, int& levels);
}

}  // namespace functionspace
}  // namespace atlas
