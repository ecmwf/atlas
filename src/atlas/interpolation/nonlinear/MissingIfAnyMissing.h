/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/interpolation/nonlinear/NonLinear.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


template <typename T>
struct MissingIfAnyMissing : NonLinear {
    bool execute( NonLinear::Matrix& W, const Field& field ) const {
        interpolation::MissingValue mv( field );
        if ( !mv ) {
            return false;
        }

        // NOTE only for scalars (for now)
        auto values        = make_view_field_values<T, 1>( field );
        auto& missingValue = mv.ref();

        // correct matrix weigths for the missing values
        // (force a missing value only if any row values is missing)
        ATLAS_ASSERT( idx_t( W.cols() ) == values.size() );

        auto data  = const_cast<Scalar*>( W.data() );
        bool modif = false;
        bool zeros = false;

        Size i = 0;
        Matrix::iterator it( W );
        for ( Size r = 0; r < W.rows(); ++r ) {
            const Matrix::iterator end = W.end( r );

            // count missing values, accumulate weights (disregarding missing values)
            size_t i_missing = i;
            size_t N_missing = 0;
            size_t N_entries = 0;

            Matrix::iterator kt( it );
            Size k = i;
            for ( ; it != end; ++it, ++i, ++N_entries ) {
                const bool miss = missingValue( values[it.col()] );

                if ( miss ) {
                    ++N_missing;
                    i_missing = i;
                }
            }

            // if any values in row are missing, force missing value
            if ( N_missing > 0 ) {
                for ( Size j = k; j < k + N_entries; ++j ) {
                    if ( j == i_missing ) {
                        data[j] = 1.;
                    }
                    else {
                        data[j] = 0.;
                        zeros   = true;
                    }
                }
                modif = true;
            }
        }

        if ( zeros && missingValue.isnan() ) {
            W.prune( 0. );
        }

        return modif;
    }

    static std::string static_type() { return "missing-if-any-missing"; }
};


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
