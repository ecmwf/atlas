/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <complex>
#include <memory>

#include "atlas/interpolation/method/sphericalvector/SparseMatrix.h"

namespace atlas {
namespace interpolation {
namespace method {
namespace detail {

using Real = double;
using Complex = std::complex<Real>;
using ComplexMatrix = SparseMatrix<Complex>;
using RealMatrix = SparseMatrix<Real>;
using Index = ComplexMatrix::Index;
using Size = ComplexMatrix::Size;
using ComplexTriplet = ComplexMatrix::Triplet;
using ComplexTriplets = ComplexMatrix::Triplets;
using RealTriplet = RealMatrix::Triplet;
using RealTriplets = RealMatrix::Triplets;
using ComplexMatPtr = std::unique_ptr<const ComplexMatrix>;
using RealMatPtr = std::unique_ptr<const RealMatrix>;


} // detail
} // method
} // interpolation
} // atlas
