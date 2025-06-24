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


#include "atlas/interpolation/nonlinear/Missing.h"

#include "eckit/types/FloatCompare.h"

#include "atlas/field/MissingValue.h"
#include "atlas/util/DataType.h"
#include "atlas/util/Metadata.h"


namespace atlas {
namespace interpolation {
namespace nonlinear {

// Factory builders
static NonLinearFactoryBuilder<MissingIfAllMissing>       __nl1(MissingIfAllMissing::static_type());
static NonLinearFactoryBuilder<MissingIfAnyMissing>       __nl2(MissingIfAnyMissing::static_type());
static NonLinearFactoryBuilder<MissingIfHeaviestMissing>  __nl3(MissingIfHeaviestMissing::static_type());

// Deprecated factory entries with "-real32" and "-real64" suffix for backwards compatibility.
// static NonLinearFactoryBuilder<MissingIfAllMissing>       __nl1_real32(MissingIfAllMissing::static_type()+"-real32");
// static NonLinearFactoryBuilder<MissingIfAnyMissing>       __nl2_real32(MissingIfAnyMissing::static_type()+"-real32");
// static NonLinearFactoryBuilder<MissingIfHeaviestMissing>  __nl3_real32(MissingIfHeaviestMissing::static_type()+"-real32");
// static NonLinearFactoryBuilder<MissingIfAllMissing>       __nl1_real64(MissingIfAllMissing::static_type()+"-real64");
// static NonLinearFactoryBuilder<MissingIfAnyMissing>       __nl2_real64(MissingIfAnyMissing::static_type()+"-real64");
// static NonLinearFactoryBuilder<MissingIfHeaviestMissing>  __nl3_real64(MissingIfHeaviestMissing::static_type()+"-real64");

namespace {
struct force_link {
    template <typename M>
    void load_builder() {
        NonLinearFactoryBuilder<M>("tmp");
    }
    force_link() {
        load_builder<MissingIfAllMissing>();
        load_builder<MissingIfAnyMissing>();
        load_builder<MissingIfHeaviestMissing>();
    }
};

std::string missing_value_type(const DataType& datatype, const field::MissingValue::Config& c) {
    std::string result;
    if (c.get("missing_value_type", result)) {
        result += "-" + datatype.str();
    }
    return result;
}

}  // namespace
void force_link_missing() {
    static force_link static_linking;
}

bool Missing::applicable(const Field& f) const {
    return field::MissingValue(f);
}

bool MissingIfAllMissing::execute(NonLinear::Matrix& W, const Field& field) const {
    return execute(W, field.array(), field.metadata());
}

bool MissingIfAllMissing::execute(NonLinear::Matrix& W, const array::Array& array, const Config& config) const {
    switch(array.datatype().kind()) {
        case (DataType::kind<double>()):        return executeT<double>(W, array, config);
        case (DataType::kind<float>()):         return executeT<float>(W, array, config);
        case (DataType::kind<int>()):           return executeT<int>(W, array, config);
        case (DataType::kind<long>()):          return executeT<long>(W, array, config);
        case (DataType::kind<unsigned long>()): return executeT<unsigned long>(W, array, config);
        default: ATLAS_NOTIMPLEMENTED;
    }
}

template<typename T>
bool MissingIfAllMissing::executeT(NonLinear::Matrix& W, const array::Array& array, const Config& config) const {
    field::MissingValue mv(missing_value_type(array.datatype(), config), config);

    auto& missingValue = mv.ref();

    ATLAS_ASSERT(array.rank() == 1);

    // NOTE only for scalars (for now)
    auto values = make_view_array_values<T, 1>(array);
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

bool MissingIfAnyMissing::execute(NonLinear::Matrix& W, const Field& field) const {
    return execute(W, field.array(), field.metadata());
}

bool MissingIfAnyMissing::execute(NonLinear::Matrix& W, const array::Array& array, const Config& config) const {
    switch(array.datatype().kind()) {
        case (DataType::kind<double>()):        return executeT<double>(W, array, config);
        case (DataType::kind<float>()):         return executeT<float>(W, array, config);
        case (DataType::kind<int>()):           return executeT<int>(W, array, config);
        case (DataType::kind<long>()):          return executeT<long>(W, array, config);
        case (DataType::kind<unsigned long>()): return executeT<unsigned long>(W, array, config);
        default: ATLAS_NOTIMPLEMENTED;
    }
}


template<typename T>
bool MissingIfAnyMissing::executeT(NonLinear::Matrix& W, const array::Array& array, const Config& config) const {
    field::MissingValue mv(missing_value_type(array.datatype(), config), config);

    auto& missingValue = mv.ref();

    // NOTE only for scalars (for now)
    auto values = make_view_array_values<T, 1>(array);
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

bool MissingIfHeaviestMissing::execute(NonLinear::Matrix& W, const Field& field) const {
    return execute(W, field.array(), field.metadata());
}

bool MissingIfHeaviestMissing::execute(NonLinear::Matrix& W, const array::Array& array, const Config& config) const {
    switch(array.datatype().kind()) {
        case (DataType::kind<double>()):        return executeT<double>(W, array, config);
        case (DataType::kind<float>()):         return executeT<float>(W, array, config);
        case (DataType::kind<int>()):           return executeT<int>(W, array, config);
        case (DataType::kind<long>()):          return executeT<long>(W, array, config);
        case (DataType::kind<unsigned long>()): return executeT<unsigned long>(W, array, config);
        default: ATLAS_NOTIMPLEMENTED;
    }
}

template<typename T>
bool MissingIfHeaviestMissing::executeT(NonLinear::Matrix& W, const array::Array& array, const Config& config) const {
    field::MissingValue mv(missing_value_type(array.datatype(), config), config);
    auto& missingValue = mv.ref();

    // NOTE only for scalars (for now)
    auto values = make_view_array_values<T, 1>(array);
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


}  // namespace nonlinear
}  // namespace interpolation
}  // namespace atlas
