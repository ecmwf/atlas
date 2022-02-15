/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

/// @file KJacobians.h
///
/// This file contains a class for working with 2D and 3D generic Jacobians.
/// The class uses eckit::maths::Matrix as a back end can be multiplied with
/// eckit::geometry::KPoint objects.

#include <iostream>

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "eckit/geometry/KPoint.h"
#include "eckit/maths/Matrix.h"

namespace atlas {
namespace util {


using BaseMatrix = eckit::maths::Matrix<double>;

template <int Rank>
class KJacobian : public BaseMatrix {
using BaseProxy = BaseMatrix::Proxy;

public:
    using BaseMatrix::Matrix;
    using BaseMatrix::operator*;

    /// @brief Default constructor.
    KJacobian() : BaseMatrix(Rank, Rank) {}

    /// @brief List constructor.
    KJacobian(std::initializer_list<std::initializer_list<double>> list) : KJacobian() {
        ATLAS_ASSERT(list.size() == Rank);
        int i = 0;
        for (const std::initializer_list<double>& subList : list) {
            ATLAS_ASSERT(subList.size() == Rank);
            int j = 0;
            for (const double& elem : subList) {
                (*this)(i, j) = elem;
                ++j;
            }
            ++i;
        }
    }

    /// @brief Base class copy constructor.
    KJacobian(const BaseMatrix& matrix) : BaseMatrix(matrix) {
#if ATLAS_BUILD_TYPE_DEBUG
        ATLAS_ASSERT(matrix.rows() == Rank);
        ATLAS_ASSERT(matrix.cols() == Rank);
#endif
    }

    /// Jacobian-KPoint multiplication.
    eckit::geometry::KPoint<Rank> operator*(const eckit::geometry::KPoint<Rank>& x) {

        auto xVec = BaseMatrix(Rank, 1);
        for (size_t i = 0; i < Rank; ++i) {
            xVec.data()[i] = x.data()[i];
        }

        const BaseMatrix yVec = this->BaseMatrix::operator*(xVec);

        return eckit::geometry::KPoint<Rank>{yVec.data()};
    }

};

/// BaseMatrix-KPoint multiplication.
template<int Rank>
eckit::geometry::KPoint<Rank> operator*(const BaseMatrix& jac, const eckit::geometry::KPoint<Rank>& x) {
#if ATLAS_BUILD_TYPE_DEBUG
        ATLAS_ASSERT(jac.rows() == Rank);
        ATLAS_ASSERT(jac.cols() == Rank);
#endif
    return KJacobian<Rank>(jac) * x;
}

using Jacobian2 = KJacobian<2>;
using Jacobian3 = KJacobian<3>;

} // namespace util
} // namespace atlas
