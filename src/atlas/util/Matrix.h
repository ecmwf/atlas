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

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "eckit/geometry/KPoint.h"
#include "eckit/maths/Matrix.h"

namespace atlas {
namespace util {


/// @brief   Matrix class.
///
/// @details Matrix class which uses eckit::maths::Matrix as a backend.
///          Matrix dimensions are strongly typed and class enables vector-
///          matrix multiplication with KPoint objects.
template <typename Value, int NRows, int NCols>
class Matrix {
public:

    using BaseType = eckit::maths::Matrix<Value, int>;

    using xPoint = eckit::geometry::KPoint<NCols>;
    using yPoint = eckit::geometry::KPoint<NRows>;

    constexpr int rows() const {
        return NRows;
    }

    constexpr int cols() const {
        return NCols;
    }

    constexpr int size() const {
        return NRows * NCols;
    }

    /// @brief Default constructor.
    Matrix() : baseMatrix_{NRows, NCols} {}

    /// @brief Data constructor (no bounds checking!).
    Matrix(const Value* data) : Matrix() {
        for (size_t i = 0; i < NRows * NCols; ++i) {
            baseMatrix_.data()[i] = data[i];
        }
    }

    /// @brief Base matrix constructor (no bounds checking in release build!).
    Matrix(const BaseType& baseMatrix) : baseMatrix_{baseMatrix} {
#if ATLAS_BUILD_TYPE_DEBUG
        ATLAS_ASSERT(baseMatrix_.rows() == NRows);
        ATLAS_ASSERT(baseMatrix_.cols() == NCols);
#endif
    }

    /// @brief List constructor (always checks bounds).
    Matrix(std::initializer_list<std::initializer_list<Value>> list) : Matrix() {
        ATLAS_ASSERT(list.size() == NRows);
        int i = 0;
        for (const std::initializer_list<Value>& subList : list) {
            ATLAS_ASSERT(subList.size() == NCols);
            int j = 0;
            for (const Value& elem : subList) {
                baseMatrix_(i, j) = elem;
                ++j;
            }
            ++i;
        }
    }

    /// @brief Get element.
    Value operator()(int i, int j) const {
        return baseMatrix_(i, j);
    }

    /// @brief Set element.
    Value& operator()(int i, int j) {
        return baseMatrix_(i, j);
    }

    /// @brief Plus operator.
    Matrix operator+(const Matrix& other) const {
        return Matrix{baseMatrix_ + other.baseMatrix()};
    }

    /// @brief Minus operator.
    Matrix operator-(const Matrix& other) const {
        return Matrix{baseMatrix_ - other.baseMatrix()};
    }

    /// @brief Multiply operator.
    template <int NCols2>
    Matrix<Value, NRows, NCols2> operator*(const Matrix<Value, NCols, NCols2>& other) const {
        return Matrix<Value, NRows, NCols2>{baseMatrix_ * other.baseMatrix()};
    }

    /// @brief Get row.
    Matrix<Value, 1, NCols> row(int i) const {
        return Matrix<Value, 1, NCols>{baseMatrix_.row(i)};
    }

    /// @brief Get col.
    Matrix<Value, NRows, 1> col(int j) const {
        return Matrix<Value, NRows, 1>{baseMatrix_.col(j)};
    }

    /// @brief Get inverse.
    Matrix inverse() const {
        return Matrix{baseMatrix_.inverse()};
    }

    /// @brief Get transpose.
    Matrix<Value, NCols, NRows> transpose() {
        // Issue: transpose() is non const in eckit/maths/MatrixLapack.h!
        return Matrix<Value, NCols, NRows>{baseMatrix_.transpose()};
    }

    /// @brief Get determinant.
    Value determinant() const {
        return baseMatrix_.determinant();
    }

    /// @brief L2 norm.
    Value norm() const {
        Value n{};
        for (size_t i = 0; i < NRows * NCols; ++i) {
            const auto elem = baseMatrix().data()[i];
            n += std::abs(elem) * std::abs(elem);
        }
        return std::sqrt(n);
    }

    /// Get signed elements of matrix (i.e., 0, +1 or -1).
    Matrix sign(Value tol = std::numeric_limits<Value>::epsilon()) const {

        const Value smallNumber = norm() * tol;
        auto sgn = Matrix{};
        for (size_t i = 0; i < NRows * NCols; ++i) {
            const auto elem = baseMatrix_.data()[i];
            sgn.baseMatrix().data()[i] = std::abs(elem) < tol ? 0. : elem < 0. ? -1. : 1.;
        }
        return sgn;
    }

    /// @brief Get const base matrix object.
    const BaseType& baseMatrix() const {
        return baseMatrix_;
    }

    /// @brief Get base matrix object.
    BaseType& baseMatrix() {
        return baseMatrix_;
    }

    /// @brief Get const matrix data pointer.
    const Value* data() const {
        return baseMatrix_.data();
    }

    /// @brief Get matrix data pointer.
    Value* data() {
        return baseMatrix_.data();
    }

    /// @brief KPoint linear transform.
    yPoint operator*(const xPoint& x) const {

        // Copy point data into column vector.
        Matrix<double, NCols, 1> xVec{x.data()};

        // Perform multiplication.
        const Matrix<double, NRows, 1> yVec = (*this) * xVec;

        // Copy column vector data to point.
        return yPoint{yVec.data()};
    }

    /// @brief Scalar mulitplication.
    Matrix operator*(Value a) const {
        auto mat = Matrix{};
        for (size_t i = 0; i < NRows * NCols; ++i) {
            mat.baseMatrix().data()[i] = baseMatrix_.data()[i] * a;
        }
        return mat;
    }

    /// @brief Equality operator.
    bool operator==(const Matrix& other) {
        for (int i = 0; i < NRows * NCols; ++i) {
            if (data()[i] != other.data()[i]) {
                return false;
            }
        }
        return true;
    }

    template <typename value, int nrows, int ncols>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<value, nrows, ncols>& mat);

private:

    BaseType baseMatrix_{};

};

template <typename Value, int NRows, int NCols>
std::ostream& operator<<(std::ostream& os, const Matrix<Value, NRows, NCols>& mat) {
    return os << mat.baseMatrix();
}

using Matrix22 = Matrix<double, 2, 2>;
using Matrix33 = Matrix<double, 3, 3>;

} // namespace util
} // namespace atlas
