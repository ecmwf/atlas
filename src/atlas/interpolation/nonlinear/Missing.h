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

#include "eckit/types/FloatCompare.h"

#include "atlas/field/MissingValue.h"
#include "atlas/interpolation/nonlinear/NonLinear.h"
#include "atlas/util/DataType.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {


struct Missing : NonLinear {
private:
    bool applicable(const Field& f) const override { return field::MissingValue(f); }
};


struct MissingIfAllMissing : Missing {
    bool execute(NonLinear::Matrix& W, const Field& field) const {
        switch(field.datatype().kind()) {
            case (DataType::kind<double>()):        return executeT<double>(W,field);
            case (DataType::kind<float>()):         return executeT<float>(W,field);
            case (DataType::kind<int>()):           return executeT<int>(W,field);
            case (DataType::kind<long>()):          return executeT<long>(W,field);
            case (DataType::kind<unsigned long>()): return executeT<unsigned long>(W,field);
            default: ATLAS_NOTIMPLEMENTED;
        }
    }
    template<typename T>
    bool executeT(NonLinear::Matrix& W, const Field& field) const {
        field::MissingValue mv(field);
        auto& missingValue = mv.ref();

        ATLAS_ASSERT(field.rank() == 1);

        auto values = make_view_field_values<T, 1>(field);
        ATLAS_ASSERT(idx_t(W.cols()) == values.shape(0));

        auto data  = const_cast<Scalar*>(W.data());
        bool modif = false;
        bool zeros = false;

        Size i = 0;
        Matrix::iterator it(W);
        for (Size r = 0; r < W.rows(); ++r) {
            const Matrix::iterator end = W.end(r);

            // count missing values, accumulate weights (disregarding missing values)
            size_t i_missing = i;
            size_t N_missing = 0;
            size_t N_entries = 0;
            Scalar sum       = 0.;

            Matrix::iterator kt(it);
            Size k = i;
            for (; it != end; ++it, ++i, ++N_entries) {
                const bool miss = missingValue(values[it.col()]);

                if (miss) {
                    ++N_missing;
                    i_missing = i;
                }
                else {
                    sum += *it;
                }
            }

            // weights redistribution: zero-weight all missing values, linear re-weighting for the others;
            // the result is missing value if all values in row are missing
            if (N_missing > 0) {
                if (N_missing == N_entries || eckit::types::is_approximately_equal(sum, 0.)) {
                    for (Size j = k; j < k + N_entries; ++j) {
                        data[j] = j == i_missing ? 1. : 0.;
                    }
                }
                else {
                    const Scalar factor = 1. / sum;
                    for (Size j = k; j < k + N_entries; ++j, ++kt) {
                        if (missingValue(values[kt.col()])) {
                            data[j] = 0.;
                            zeros   = true;
                        }
                        else {
                            data[j] *= factor;
                        }
                    }
                }
                modif = true;
            }
        }

        if (zeros && missingValue.isnan()) {
            W.prune(0.);
        }

        return modif;
    }

    static std::string static_type() { return "missing-if-all-missing"; }
};


struct MissingIfAnyMissing : Missing {
    bool execute(NonLinear::Matrix& W, const Field& field) const {
        switch(field.datatype().kind()) {
            case (DataType::kind<double>()):        return executeT<double>(W,field);
            case (DataType::kind<float>()):         return executeT<float>(W,field);
            case (DataType::kind<int>()):           return executeT<int>(W,field);
            case (DataType::kind<long>()):          return executeT<long>(W,field);
            case (DataType::kind<unsigned long>()): return executeT<unsigned long>(W,field);
            default: ATLAS_NOTIMPLEMENTED;
        }
    }

    template<typename T>
    bool executeT(NonLinear::Matrix& W, const Field& field) const {
        field::MissingValue mv(field);
        auto& missingValue = mv.ref();

        // NOTE only for scalars (for now)
        auto values = make_view_field_values<T, 1>(field);
        ATLAS_ASSERT(idx_t(W.cols()) == values.size());

        auto data  = const_cast<Scalar*>(W.data());
        bool modif = false;
        bool zeros = false;

        Size i = 0;
        Matrix::iterator it(W);
        for (Size r = 0; r < W.rows(); ++r) {
            const Matrix::iterator end = W.end(r);

            // count missing values, accumulate weights (disregarding missing values)
            size_t i_missing = i;
            size_t N_missing = 0;
            size_t N_entries = 0;

            Matrix::iterator kt(it);
            Size k = i;
            for (; it != end; ++it, ++i, ++N_entries) {
                const bool miss = missingValue(values[it.col()]);

                if (miss) {
                    ++N_missing;
                    i_missing = i;
                }
            }

            // if any values in row are missing, force missing value
            if (N_missing > 0) {
                for (Size j = k; j < k + N_entries; ++j) {
                    if (j == i_missing) {
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

        if (zeros && missingValue.isnan()) {
            W.prune(0.);
        }

        return modif;
    }

    static std::string static_type() { return "missing-if-any-missing"; }
};


struct MissingIfHeaviestMissing : Missing {
    bool execute(NonLinear::Matrix& W, const Field& field) const {
        switch(field.datatype().kind()) {
            case (DataType::kind<double>()):        return executeT<double>(W,field);
            case (DataType::kind<float>()):         return executeT<float>(W,field);
            case (DataType::kind<int>()):           return executeT<int>(W,field);
            case (DataType::kind<long>()):          return executeT<long>(W,field);
            case (DataType::kind<unsigned long>()): return executeT<unsigned long>(W,field);
            default: ATLAS_NOTIMPLEMENTED;
        }
    }

    template<typename T>
    bool executeT(NonLinear::Matrix& W, const Field& field) const {
        field::MissingValue mv(field);
        auto& missingValue = mv.ref();

        // NOTE only for scalars (for now)
        auto values = make_view_field_values<T, 1>(field);
        ATLAS_ASSERT(idx_t(W.cols()) == values.size());

        auto data  = const_cast<Scalar*>(W.data());
        bool modif = false;
        bool zeros = false;

        Size i = 0;
        Matrix::iterator it(W);
        for (Size r = 0; r < W.rows(); ++r) {
            const Matrix::iterator end = W.end(r);

            // count missing values, accumulate weights (disregarding missing values) and find maximum weight in row
            size_t i_missing         = i;
            size_t N_missing         = 0;
            size_t N_entries         = 0;
            Scalar sum               = 0.;
            Scalar heaviest          = -1.;
            bool heaviest_is_missing = false;

            Matrix::iterator kt(it);
            Size k = i;
            for (; it != end; ++it, ++i, ++N_entries) {
                const bool miss = missingValue(values[it.col()]);

                if (miss) {
                    ++N_missing;
                    i_missing = i;
                }
                else {
                    sum += *it;
                }

                if (heaviest < data[i]) {
                    heaviest            = data[i];
                    heaviest_is_missing = miss;
                }
            }

            // weights redistribution: zero-weight all missing values, linear re-weighting for the others;
            // if all values are missing, or the closest value is missing, force missing value
            if (N_missing > 0) {
                if (N_missing == N_entries || heaviest_is_missing || eckit::types::is_approximately_equal(sum, 0.)) {
                    for (Size j = k; j < k + N_entries; ++j) {
                        data[j] = j == i_missing ? 1. : 0.;
                    }
                }
                else {
                    const Scalar factor = 1. / sum;
                    for (Size j = k; j < k + N_entries; ++j, ++kt) {
                        if (missingValue(values[kt.col()])) {
                            data[j] = 0.;
                            zeros   = true;
                        }
                        else {
                            data[j] *= factor;
                        }
                    }
                }
                modif = true;
            }
        }

        if (zeros && missingValue.isnan()) {
            W.prune(0.);
        }

        return modif;
    }

    static std::string static_type() { return "missing-if-heaviest-missing"; }
};

}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
