/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

/// @file KJacobian.h
///
/// Contains functions for performing Jacobian-KPoint multiplications.
/// Uses eckit::maths::Matrix type to represent Jacobian.


#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "eckit/geometry/KPoint.h"
#include "eckit/maths/Matrix.h"

namespace atlas {
namespace util {

using KJacobian = eckit::maths::Matrix<double>;

/// @brief Create a Jacobian from an initialiser list.
/// @note  To be removed once constructor is added to eckit::maths::Matrix.
template<int Rank>
KJacobian make_KJacobian(std::initializer_list<std::initializer_list<double>> list) {
    auto jac = KJacobian(Rank, Rank);
    ATLAS_ASSERT(list.size() == Rank);
    int i = 0;
    for (const std::initializer_list<double>& subList : list) {
        ATLAS_ASSERT(subList.size() == Rank);
        int j = 0;
        for (const double& elem : subList) {
            jac(i, j) = elem;
            ++j;
        }
        ++i;
    }
    return jac;
}

/// Jacobian-KPoint multiplication.
template<int Rank>
eckit::geometry::KPoint<Rank> operator*(const KJacobian& jac, const eckit::geometry::KPoint<Rank>& x) {
#if ATLAS_BUILD_TYPE_DEBUG
        ATLAS_ASSERT(jac.rows() == Rank);
        ATLAS_ASSERT(jac.cols() == Rank);
#endif

        auto xVec = eckit::maths::ColVector<double>(Rank);
        std::copy(x.data(), x.data() + Rank, xVec.data());

        const eckit::maths::ColVector<double> yVec = jac * xVec;

        return eckit::geometry::KPoint<Rank>{yVec.data()};
}

const auto make_Jacobian2 = make_KJacobian<2>;
const auto make_Jacobian3 = make_KJacobian<3>;

} // namespace util
} // namespace atlas
