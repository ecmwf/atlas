/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "SpectralInterface.h"

#include "atlas/array/LocalView.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/runtime/Exception.h"
#include "atlas/trans/Trans.h"

namespace atlas {
namespace functionspace {

namespace detail {
struct SpectralFortranAccess {
    const Spectral& fs_;
    SpectralFortranAccess(const Spectral& fs): fs_(fs) {}
    int nump() const { return fs_.nump(); }
    array::LocalView<const int, 1> nvalue() const { return fs_.nvalue(); }
    array::LocalView<const int, 1> nmyms() const { return fs_.nmyms(); }
    array::LocalView<const int, 1> nasm0() const { return fs_.nasm0(); }
};
}  // namespace detail

using detail::SpectralFortranAccess;

// ----------------------------------------------------------------------

extern "C" {
const detail::Spectral* atlas__SpectralFunctionSpace__new__config(const eckit::Configuration* config) {
    ATLAS_ASSERT(config != nullptr);
    return new detail::Spectral(*config);
}

void atlas__SpectralFunctionSpace__delete(detail::Spectral* This) {
    ATLAS_ASSERT(This != nullptr);
    delete This;
}

field::FieldImpl* atlas__fs__Spectral__create_field(const detail::Spectral* This, const eckit::Configuration* options) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(options);
    field::FieldImpl* field;
    {
        Field f = This->createField(*options);
        field   = f.get();
        field->attach();
    }
    field->detach();
    return field;
}

void atlas__SpectralFunctionSpace__gather(const detail::Spectral* This, const field::FieldImpl* local,
                                          field::FieldImpl* global) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(global != nullptr);
    ATLAS_ASSERT(local != nullptr);
    const Field l(local);
    Field g(global);
    This->gather(l, g);
}

void atlas__SpectralFunctionSpace__scatter(const detail::Spectral* This, const field::FieldImpl* global,
                                           field::FieldImpl* local) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(global != nullptr);
    ATLAS_ASSERT(local != nullptr);
    const Field g(global);
    Field l(local);
    This->scatter(g, l);
}

void atlas__SpectralFunctionSpace__gather_fieldset(const detail::Spectral* This, const field::FieldSetImpl* local,
                                                   field::FieldSetImpl* global) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(global != nullptr);
    ATLAS_ASSERT(local != nullptr);
    const FieldSet l(local);
    FieldSet g(global);
    This->gather(l, g);
}

void atlas__SpectralFunctionSpace__scatter_fieldset(const detail::Spectral* This, const field::FieldSetImpl* global,
                                                    field::FieldSetImpl* local) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(global != nullptr);
    ATLAS_ASSERT(local != nullptr);
    const FieldSet g(global);
    FieldSet l(local);
    This->scatter(g, l);
}

void atlas__SpectralFunctionSpace__norm(const detail::Spectral* This, const field::FieldImpl* field, double norm[],
                                        int rank) {
    ATLAS_ASSERT(This != nullptr);
    ATLAS_ASSERT(field != nullptr);
    ATLAS_ASSERT(norm != nullptr);
    This->norm(field, norm, rank);
}

void atlas__SpectralFunctionSpace__nspec2(const detail::Spectral* This, int& nspec2) {
    nspec2 = This->nb_spectral_coefficients();
}

void atlas__SpectralFunctionSpace__nspec2g(const detail::Spectral* This, int& nspec2g) {
    nspec2g = This->nb_spectral_coefficients_global();
}

void atlas__SpectralFunctionSpace__truncation(const detail::Spectral* This, int& truncation) {
    truncation = This->truncation();
}

void atlas__SpectralFunctionSpace__nump(const detail::Spectral* This, int& nump) {
    nump = detail::SpectralFortranAccess(*This).nump();
}

void atlas__SpectralFunctionSpace__nmyms(const detail::Spectral* This, const int*& nmyms, int& size) {
    const auto nmyms_ = SpectralFortranAccess(*This).nmyms();
    nmyms             = nmyms_.data();
    size              = nmyms_.size();
}

void atlas__SpectralFunctionSpace__nasm0(const detail::Spectral* This, const int*& nasm0, int& size) {
    const auto nasm0_ = SpectralFortranAccess(*This).nasm0();
    nasm0             = nasm0_.data();
    size              = nasm0_.size();
}

void atlas__SpectralFunctionSpace__nvalue(const detail::Spectral* This, const int*& nvalue, int& size) {
    const auto nvalue_ = SpectralFortranAccess(*This).nvalue();
    nvalue             = nvalue_.data();
    size               = nvalue_.size();
}

void atlas__SpectralFunctionSpace__levels(const detail::Spectral* This, int& levels) {
    levels = This->levels();
}


}  // extern "C"

}  // namespace functionspace
}  // namespace atlas
