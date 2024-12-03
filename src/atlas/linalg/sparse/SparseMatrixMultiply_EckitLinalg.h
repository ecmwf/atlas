/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */


#pragma once

#include "atlas/linalg/sparse/SparseMatrixMultiply.h"

namespace atlas {
namespace linalg {
namespace sparse {

template <>
struct SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 1, double, const double, double> {
    static void multiply(const SparseMatrixView<double>&, const View<const double, 1>& src, View<double, 1>& tgt,
                         const Configuration&);
    static void multiply_add(const SparseMatrixView<double>&, const View<const double, 1>& src, View<double, 1>& tgt,
                             const Configuration&);
};

template <>
struct SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_right, 2, double, const double, double> {
    static void multiply(const SparseMatrixView<double>&, const View<const double, 2>& src, View<double, 2>& tgt,
                         const Configuration&);
    static void multiply_add(const SparseMatrixView<double>&, const View<const double, 2>& src, View<double, 2>& tgt,
                             const Configuration&);
};


template <>
struct SparseMatrixMultiply<backend::eckit_linalg, Indexing::layout_left, 1, double, const double, double> {
    static void multiply(const SparseMatrixView<double>&, const View<const double, 1>& src, View<double, 1>& tgt,
                         const Configuration&);
    static void multiply_add(const SparseMatrixView<double>&, const View<const double, 1>& src, View<double, 1>& tgt,
                             const Configuration&);
};

}  // namespace sparse
}  // namespace linalg
}  // namespace atlas
