/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

/// @file Matrix.h
///
/// This file contains classes and functions for working with small square matrices.

#include <iostream>

#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "eckit/geometry/KPoint.h"
#include "eckit/maths/Matrix.h"

#include "atlas/util/Point.h"

namespace atlas {
namespace util {


/// @brief   Matrix class.
///
/// @details Matrix class which uses eckit::maths::Matrix as a backend.
///          Matrix dimensions are strongly typed and class enables vector-
///          matrix multiplication with KPoint classes.
template <typename Value, int NRows, int NCols>
class Matrix {
public:

    using BaseType = eckit::maths::Matrix<Value, int>;

    using xPoint = eckit::geometry::KPoint<NRows>;
    using yPoint = eckit::geometry::KPoint<NCols>;


    /// @brief Default constructor.
    Matrix() : baseMatrix_{NRows, NCols} {}

    /// @brief base matrix constructor.
    Matrix(const BaseType& baseMatrix) : baseMatrix_{baseMatrix} {
#if ATLAS_BUILD_TYPE_DEBUG
        ATLAS_ASSERT(baseMatrix_.rows() == NRows);
        ATLAS_ASSERT(baseMatrix_.cols() == NCols);
#endif
    }

    /// @brief List constructor.
    Matrix(std::initializer_list<std::initializer_list<Value>> list) : Matrix() {
        // Get pointer to first element of data_.
        int i = 0;
        for (const std::initializer_list<Value>& subList : list) {
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
    Matrix transpose() const {
        return Matrix{baseMatrix_.transpose()};
    }

    /// @brief Get determinant.
    Value determinant() const {
        return baseMatrix_.determinant();
    }

    /// @brief L2 norm.
    Value norm() const {
        Value n{};
        for (size_t i = 0; i < NRows * NCols; ++i) {
            const auto elem = baseMatrix_.data()[i];
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
            sgn.baseMatrix_.data()[i] = std::abs(elem) < tol ? 0. : elem < 0. ? -1. : 1.;
        }
        return sgn;
    }

    /// @brief Get base matrix object.
    const BaseType& baseMatrix() const {
        return baseMatrix_;
    }

    /// @brief Get matrix data pointer.
    const Value* data() const {
        return baseMatrix().data();
    }

    /// @brief KPoint linear transform.
    yPoint operator*(const xPoint& x) const {

        // Const cast needed due to issues with ConstProxy type.
        auto xPtr = const_cast<double*>(x.data());
        const Matrix<double, NRows, 1> xVec{typename BaseType::Proxy(xPtr, NRows, 1)};
        const Matrix<double, NCols, 1> yVec = (*this) * xVec;

        return yPoint{yVec.data()};
    }

    /// @brief Scalar mulitplication.
    Matrix operator*(Value a) const {
        auto mat = Matrix{};
        for (size_t i = 0; i < NRows * NCols; ++i) {
            mat.baseMatrix_.data()[i] = baseMatrix_.data()[i] * a;
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
    friend std::ostream& operator <<(std::ostream& os, const Matrix<value, nrows, ncols>& mat);

private:

    BaseType baseMatrix_{};

};

template <typename Value, int NRows, int NCols>
std::ostream& operator <<(std::ostream& os, const Matrix<Value, NRows, NCols>& mat) {
    return os << mat.baseMatrix_;
}

using Matrix22 = Matrix<double, 2, 2>;
using Matrix33 = Matrix<double, 3, 3>;

} // namespace util




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
