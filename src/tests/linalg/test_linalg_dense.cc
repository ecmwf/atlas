/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <memory>
#include <tuple>
#include <vector>

#include "atlas/library/config.h"
#if ATLAS_ECKIT_HAVE_ECKIT_585
#include "eckit/linalg/LinearAlgebraDense.h"
#else
#include "eckit/linalg/LinearAlgebra.h"
#endif

#include "eckit/linalg/Matrix.h"
#include "eckit/linalg/Vector.h"

#include "atlas/linalg/dense.h"
#include "atlas/linalg/View.h"
#include "atlas/array/ArrayView.h"
#include "atlas/array/Array.h"
#include "atlas/array/MakeView.h"

#include "tests/AtlasTestEnvironment.h"


using namespace atlas::linalg;

namespace atlas {
namespace test {

//----------------------------------------------------------------------------------------------------------------------

// strings to be used in the tests
static std::string eckit_linalg = dense::backend::eckit_linalg::type();

//----------------------------------------------------------------------------------------------------------------------

// Only reason to define these derived classes is for nicer constructors and convenience in the tests

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

CASE("test_linalg_make_view_readonly_elements") {
    // Linalg conversion must retain pointee constness even for mutable view objects.
    using ReadonlyArray = array::ArrayView<const double, 1>;
    using ReadonlyLocal = array::LocalView<const double, 1>;
    static_assert(std::is_same_v<typename linalg::view_type<ReadonlyArray>::element_type, const double>);
    static_assert(std::is_same_v<typename linalg::view_type<ReadonlyLocal>::element_type, const double>);

    std::unique_ptr<array::Array> storage{array::Array::create<double>(3)};
    auto writable_array = array::make_host_view<double, 1>(*storage);
    writable_array(0) = 1.;
    writable_array(1) = 2.;
    writable_array(2) = 3.;
    auto* data = writable_array.data();
    array::ArrayShape shape{3};
    array::LocalView<double, 1> writable_source{data, shape};
    auto readonly_local = writable_source.as_const();
    auto local_result = linalg::make_view(readonly_local);
    static_assert(!std::is_assignable_v<decltype(local_result(1)), double>);
    EXPECT(local_result.data() == data);
    EXPECT(local_result.shape(0) == 3);
    EXPECT(local_result(1) == 2.);

    auto readonly_array = array::make_host_view<const double, 1>(*storage);
    auto array_result = linalg::make_view(readonly_array);
    static_assert(!std::is_assignable_v<decltype(array_result(1)), double>);
    EXPECT(array_result.data() == data);
    EXPECT(array_result(2) == 3.);

    auto writable_result = linalg::make_view(writable_source);
    static_assert(std::is_same_v<typename decltype(writable_result)::element_type, double>);
    writable_result(1) = 42.;
    EXPECT(local_result(1) == 42.);
    EXPECT(array_result(1) == 42.);

    auto const_result = linalg::make_view(std::as_const(writable_source));
    static_assert(std::is_same_v<typename decltype(const_result)::element_type, const double>);
    EXPECT(const_result(1) == 42.);
}

CASE("test configuration via resource") {
    if (atlas::Library::instance().linalgDenseBackend().empty()) {
#if ATLAS_ECKIT_HAVE_ECKIT_585
        if (eckit::linalg::LinearAlgebraDense::hasBackend("mkl")) {
#else
        if (eckit::linalg::LinearAlgebra::hasBackend("mkl")) {
#endif
            EXPECT_EQ(dense::current_backend().type(), "mkl");
        }
        else {
            EXPECT_EQ(dense::current_backend().type(), eckit_linalg);
        }
    }
    else {
        EXPECT_EQ(dense::current_backend().type(), atlas::Library::instance().linalgDenseBackend());
    }
}

CASE("test backend functionalities") {
    EXPECT_EQ(dense::current_backend().getString("backend", "undefined"), "undefined");

    dense::current_backend(eckit_linalg);
    EXPECT_EQ(dense::current_backend().type(), "eckit_linalg");
    EXPECT_EQ(dense::current_backend().getString("backend", "undefined"), "undefined");
    dense::current_backend().set("backend", "default");
    EXPECT_EQ(dense::current_backend().getString("backend"), "default");

    EXPECT_EQ(dense::default_backend(eckit_linalg).getString("backend"), "default");

    dense::default_backend(eckit_linalg).set("backend", "generic");
    EXPECT_EQ(dense::default_backend(eckit_linalg).getString("backend"), "generic");

    const dense::Backend backend_default      = dense::Backend();
    const dense::Backend backend_eckit_linalg = dense::backend::eckit_linalg();
    EXPECT_EQ(backend_default.type(), eckit_linalg);
    EXPECT_EQ(backend_eckit_linalg.type(), eckit_linalg);

    EXPECT_EQ(std::string(backend_eckit_linalg), eckit_linalg);
}

//----------------------------------------------------------------------------------------------------------------------

CASE("matrix matrix multiply (gemm)") {
    Matrix A = Matrix{{1., -2.}, {-4., 2.}};
    Matrix Z = Matrix{{9., -6.}, {-12., 12.}};

    SECTION("bad sizes") {
        Matrix Y = Matrix{{-42., -42.}, {-42., -42.}};
        EXPECT_THROWS_AS(linalg::matrix_multiply(A, Matrix{{0, 0}}, Y), eckit::AssertionFailed);
    }

    std::vector<std::string> backends{"eckit_linalg", "generic", "openmp", "lapack", "eigen"};
    for (auto& backend : backends) {
        if (dense::Backend{backend}.available()) {
            SECTION(backend) {
                Matrix Y = Matrix{{-42., -42.}, {-42., -42.}};
                linalg::matrix_multiply(A, A, Y, dense::Backend{backend});
                expect_equal(Y, Z);
            }
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------

}  // namespace test
}  // namespace atlas

int main(int argc, char** argv) {
    return atlas::test::run(argc, argv);
}
