/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <vector>

#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/array.h"
#include "atlas/linalg/sparse.h"

using namespace atlas::linalg;

namespace atlas {
namespace test {

class Vector : public eckit::linalg::Vector {
public:
    using Scalar = eckit::linalg::Scalar;
    using eckit::linalg::Vector::Vector;
    Vector(const std::initializer_list<Scalar>& v): eckit::linalg::Vector::Vector(v.size()) {
        size_t i = 0;
        for (auto& s : v) {
            operator[](i++) = s;
        }
    }
};

class Matrix : public eckit::linalg::Matrix {
public:
    using Scalar = eckit::linalg::Scalar;
    using eckit::linalg::Matrix::Matrix;

    Matrix(const std::initializer_list<std::vector<Scalar>>& m):
        eckit::linalg::Matrix::Matrix(m.size(), m.size() ? m.begin()->size() : 0) {
        size_t r = 0;
        for (auto& row : m) {
            for (size_t c = 0; c < cols(); ++c) {
                operator()(r, c) = row[c];
            }
            ++r;
        }
    }
};

// 2D array constructable from eckit::linalg::Matrix
// Indexing/memorylayout and data type can be customized for testing
template <typename Value, Indexing indexing = Indexing::layout_left>
struct ArrayMatrix {
    array::ArrayView<Value, 2> view() {
        return host_view();
    }
    array::ArrayView<Value, 2> host_view() {
        array.updateHost();
        return array::make_host_view<Value, 2>(array);
    }
    array::ArrayView<Value, 2> device_view() {
        array.updateDevice();
        return array::make_device_view<Value, 2>(array);
    }
    ArrayMatrix(const eckit::linalg::Matrix& m): ArrayMatrix(m.rows(), m.cols()) {
        auto view_ = array::make_view<Value, 2>(array);
        for (int r = 0; r < m.rows(); ++r) {
            for (int c = 0; c < m.cols(); ++c) {
                auto& v = layout_left ? view_(r, c) : view_(c, r);
                v       = m(r, c);
            }
        }
    }
    ArrayMatrix(int r, int c): array(make_shape(r, c)) {}

private:
    static constexpr bool layout_left = (indexing == Indexing::layout_left);
    static array::ArrayShape make_shape(int rows, int cols) {
        return layout_left ? array::make_shape(rows, cols) : array::make_shape(cols, rows);
    }
    array::ArrayT<Value> array;
};

// 1D array constructable from eckit::linalg::Vector
template <typename Value>
struct ArrayVector {
    array::ArrayView<Value, 1> view() {
        return host_view();
    }
    array::ArrayView<Value, 1> host_view() {
        array.updateHost();
        return array::make_host_view<Value, 1>(array);
    }
    array::ArrayView<Value, 1> device_view() {
        array.updateDevice();
        return array::make_device_view<Value, 1>(array);
    }
    ArrayVector(const eckit::linalg::Vector& v): ArrayVector(v.size()) {
        auto view_ = array::make_view<Value, 1>(array);
        for (int n = 0; n < v.size(); ++n) {
            view_[n] = v[n];
        }
    }
    ArrayVector(int size): array(size) {}

private:
    array::ArrayT<Value> array;
};

//----------------------------------------------------------------------------------------------------------------------

template <typename T>
void expect_equal(T* v, T* r, size_t s) {
    EXPECT(is_approximately_equal(eckit::testing::make_view(v, s), eckit::testing::make_view(r, s), T(1.e-5)));
}

template <class T1, class T2>
void expect_equal(const T1& v, const T2& r) {
    expect_equal(v.data(), r.data(), std::min(v.size(), r.size()));
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas
