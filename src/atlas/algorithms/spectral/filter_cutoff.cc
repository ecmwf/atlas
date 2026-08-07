/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/algorithms/spectral/filter_cutoff.h"

#include "atlas/array/MakeView.h"
#include "atlas/field/detail/FieldImpl.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/runtime/Exception.h"

namespace atlas {
namespace spectral {

void filter_cutoff(const functionspace::Spectral& spectral, Field& field, int cutoff) {
    if (cutoff < 0) {
        throw_Exception("Spectral cutoff must be non-negative", Here());
    }

    auto implementation_for_datatype = [&](auto datatype) {
        using Value = decltype(datatype);
        constexpr Value zero{0};
        if (field.levels() == 0) {
            auto coefficients = array::make_view<Value, 1>(field);
            spectral.parallel_for([&](idx_t real, idx_t imag, int n, int) {
                if (n > cutoff) {
                    coefficients(real) = zero;
                    coefficients(imag) = zero;
                }
            });
        }
        else {
            auto coefficients = array::make_view<Value, 2>(field);
            const int levels = field.levels();
            spectral.parallel_for([&](idx_t real, idx_t imag, int n, int) {
                if (n > cutoff) {
                    for (int level = 0; level < levels; ++level) {
                        coefficients(real, level) = zero;
                        coefficients(imag, level) = zero;
                    }
                }
            });
        }
    };

    switch (field.datatype().kind()) {
        case array::DataType::KIND_REAL32: implementation_for_datatype(float{}); break;
        case array::DataType::KIND_REAL64: implementation_for_datatype(double{}); break;
        default:
            throw_Exception("Spectral cutoff only supports real32 and real64 spectral fields", Here());
    }

    field.set_dirty();
}

}  // namespace spectral
}  // namespace atlas

extern "C" {
void atlas__spectral__filter_cutoff(const atlas::functionspace::detail::Spectral* spectral, atlas::field::FieldImpl* field,
                                    int cutoff) {
    ATLAS_ASSERT(spectral != nullptr);
    ATLAS_ASSERT(field != nullptr);

    atlas::Field atlas_field(field);
    atlas::spectral::filter_cutoff(atlas::functionspace::Spectral(atlas::FunctionSpace(spectral)), atlas_field, cutoff);
}
}
