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
class Matrix2 {
public:

    /// @brief Default constructor.
    Matrix2() = default;

    /// @brief Data constructor.
    Matrix2(const std::array<Value, 4>& a) : data_{a} {}

    /// @brief Elements constructor.
    Matrix2(Value a00, Value a01, Value a10, Value a11) : data_{{a00, a01, a10, a11}} {}

    /// @brief return data array.
    std::array<Value, 4>& data() {return data_;}

    /// @brief return const data array.
    const std::array<Value, 4>& data() const {return data_;}

    /// @brief get element.
    const Value operator()(int i, int j) const {return data_[2 * i + j];}

    /// @brief set element.
    Value& operator()(int i, int j) {return data_[2 * i + j];}

    /// @brief Determinant of matrix.
    Value det() const {
        return data_[0] * data_[3] - data_[1] * data_[2];
    }

    /// @brief matrix-scalar multiplication.
    Matrix2 operator*(Value a) const {
        return Matrix2{data_[0] * a, data_[1] * a, data_[2] * a, data_[3] * a};
    }

    /// @brief matrix-vector multiplication.
    Point2 operator*(const Point2& x) const {
        return Point2{x[0] * data_[0] + data_[1] * x[1], x[0] * data_[2] + data_[3] * x[1]};
    }

    /// @brief matrix-matrix multiplication.
    Matrix2 operator*(const Matrix2<Value>& B) const {
        return Matrix2<Value>{data_[0] * B.data_[0] + data_[1] * B.data_[2],
                              data_[0] * B.data_[1] + data_[1] * B.data_[3],
                              data_[2] * B.data_[0] + data_[3] * B.data_[2],
                              data_[2] * B.data_[1] + data_[3] * B.data_[3]};
    }

    /// @brief Inverse matrix.
    Matrix2 inverse() const {
        return Matrix2<Value>{data_[3], -data_[1], -data_[2], data_[0]} * (1. / det());
    }

    /// @brief Get signed elements of matrix (i.e., 0, +1 or -1).
    Matrix2 sign() const {
        const Value smallNumber = det() * std::numeric_limits<Value>::epsilon();
        const auto signValue     = [&](Value number) -> Value {
            return std::abs(number) < smallNumber ? 0. : number < 0. ? -1. : 1.;
        };
        return Matrix2<Value>{signValue(data_[0]), signValue(data_[1]), signValue(data_[2]), signValue(data_[3])};
    }

private:
    // Data storage.
    std::array<Value, 4> data_{};

};

/// @brief ostream insertion operator.
template <typename Value>
inline std::ostream& operator<<(std::ostream& out, const Matrix2<Value>& A) {

    out << std::to_string(A(0, 0)) + " " +
           std::to_string(A(0, 1)) + "\n" +
           std::to_string(A(1, 0)) + " " +
           std::to_string(A(1, 1)) + "\n";

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
class JacobianXY : public Matrix2<double> {
public:
using Matrix2<double>::Matrix2;


    /// @brief Converting copy constructor.
    JacobianXY(const Matrix2<double>& mat) : Matrix2<double>(mat) {};

    /// @brief Converting move constructor.
    JacobianXY(Matrix2<double>&& mat) : Matrix2<double>(mat) {};

    /// @brief Converting copy assignment.
    JacobianXY& operator=(const Matrix2<double>& mat) {
        this->data() = mat.data();
        return *this;
    }

    /// @brief Converting move assignment.
    JacobianXY& operator=(Matrix2<double>&& mat) {
        this->data() = std::move(mat.data());
        return *this;
    }

    /// @brief Finite difference constructor.
    ///
    /// @param f00 = f(X0      , X1      )
    /// @param f10 = f(X0 + dx0, X1      )
    /// @param f10 = f(X0,     , X1 + dx1)
    JacobianXY(const Point2& f00, const Point2& f10, const Point2& f01, double dx0 = 1., double dx1 = 1.) :
        Matrix2{(f10[0] - f00[0]) / dx0, (f01[0] - f00[0]) / dx1,
                (f10[1] - f00[1]) / dx0, (f01[1] - f00[1]) / dx1} {}
};


} // namespace atlas
