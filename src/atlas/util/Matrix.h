/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

/// @file Matrix.h
///
/// This file contains classes and functions for working with small matrices.

#include <iostream>

#include "atlas/util/Point.h"

namespace atlas {

/// @brief   Simple 2x2 matrix with static memory allocation.
///
/// @details Matrix class which can be used in conjunction with other matrices
///          and Point2 classes.

template <typename Value>
class SquareMatrix2 {
    // Data storage.
    std::array<std::array<Value, 2>, 2> data_{};

public:

    /// @brief Default constructor.
    SquareMatrix2() = default;

    /// @brief List constructor.
    SquareMatrix2(std::initializer_list<std::initializer_list<Value>> list)  {
        // Get pointer to first element of data_.
        Value* elemPtr = data_.front().data();
        for (const std::initializer_list<Value>& subList : list) {
            for (const Value& elem : subList) {
                *elemPtr++ = elem;
            }
        }
    }

    /// @brief std::array constructor.
    SquareMatrix2(const std::array<std::array<Value, 2>, 2>& arr) : data_{arr} {}

    /// @brief return data array.
    std::array<std::array<Value, 2>, 2>& data() {return data_;}

    /// @brief return const data array.
    const std::array<std::array<Value, 2>, 2> data() const {return data_;}

    /// @brief get row.
    const std::array<Value, 2>& operator[](size_t i) const {return data_[i];}

    /// @brief set row.
    std::array<Value, 2>& operator[](size_t i) {return data_[i];}

    /// @brief Determinant of matrix.
    Value det() const {
        return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
    }

    /// @brief matrix-scalar multiplication.
    SquareMatrix2 operator*(Value a) const {
        return SquareMatrix2{{data_[0][0] * a, data_[0][1] * a},
                             {data_[1][0] * a, data_[1][1] * a}};
    }

    /// @brief matrix-vector multiplication.
    Point2 operator*(const Point2& x) const {
        return Point2{data_[0][0] * x[0] + data_[0][1] * x[1], data_[1][0] * x[0] + data_[1][1] * x[1]};
    }

    /// @brief matrix-matrix multiplication.
    SquareMatrix2 operator*(const SquareMatrix2<Value>& B) const {
        return SquareMatrix2<Value>{{data_[0][0] * B[0][0] + data_[0][1] * B[1][0],
                                     data_[0][0] * B[0][1] + data_[0][1] * B[1][1]},
                                    {data_[1][0] * B[0][0] + data_[1][1] * B[1][0],
                                     data_[1][0] * B[0][1] + data_[1][1] * B[1][1]}};
    }

    /// @brief Inverse matrix.
    SquareMatrix2 inverse() const {
        return SquareMatrix2<Value>{{data_[1][1], -data_[0][1]},
                                    {-data_[1][0], data_[0][0]}} * (1. / det());
    }

    /// @brief Get signed elements of matrix (i.e., 0, +1 or -1).
    SquareMatrix2 sign() const {
        const Value smallNumber = det() * std::numeric_limits<Value>::epsilon();
        const auto signValue     = [&](Value number) -> Value {
            return std::abs(number) < smallNumber ? 0. : number < 0. ? -1. : 1.;
        };
        return SquareMatrix2<Value>{{signValue(data_[0][0]), signValue(data_[0][1])},
                                    {signValue(data_[1][0]), signValue(data_[1][1])}};
    }



};

/// @brief ostream insertion operator.
template <typename Value>
inline std::ostream& operator<<(std::ostream& out, const SquareMatrix2<Value>& A) {

    out << "{{" << A[0][0] << "," << A[0][1] << "},{"
                << A[1][0] << "," << A[1][1] << "}}";

    return out;
}


/// @brief   Jacobian class for a 2D vector field.
///
/// @details Includes a finite difference constructor which takes three f(x)
///          values arranged as follows:
///
///             ^
///             | *f(X0, X1 + dx1)
///             |
///          x1 |
///             |
///             | *f(X0, X1)  *f(X0 + dx0, X1)
///             +---------------------------->
///                          x0
class JacobianXY : public SquareMatrix2<double> {
public:
using SquareMatrix2<double>::SquareMatrix2;


    /// @brief Converting copy constructor.
    JacobianXY(const SquareMatrix2<double>& mat) : SquareMatrix2<double>(mat) {};

    /// @brief Converting move constructor.
    JacobianXY(SquareMatrix2<double>&& mat) : SquareMatrix2<double>(mat) {};

    /// @brief Converting copy assignment.
    JacobianXY& operator=(const SquareMatrix2<double>& mat) {
        this->data() = mat.data();
        return *this;
    }

    /// @brief Converting move assignment.
    JacobianXY& operator=(SquareMatrix2<double>&& mat) {
        this->data() = std::move(mat.data());
        return *this;
    }

    /// @brief Finite difference constructor.
    ///
    /// @param f00 = f(X0      , X1      )
    /// @param f10 = f(X0 + dx0, X1      )
    /// @param f01 = f(X0,     , X1 + dx1)
    JacobianXY(const Point2& f00, const Point2& f10, const Point2& f01, double dx0 = 1., double dx1 = 1.) :
        SquareMatrix2{{(f10[0] - f00[0]) / dx0, (f01[0] - f00[0]) / dx1},
                      {(f10[1] - f00[1]) / dx0, (f01[1] - f00[1]) / dx1}} {}
};

} // namespace atlas
