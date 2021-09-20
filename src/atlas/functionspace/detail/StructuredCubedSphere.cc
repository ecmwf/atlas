/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/functionspace/detail/StructuredCubedSphere.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace functionspace {
namespace detail {

namespace  {
#if ATLAS_BUILD_TYPE_DEBUG
  constexpr bool checkBounds = true;
#else
  constexpr bool checkBounds = false;
#endif
}

StructuredCubedSphere::~StructuredCubedSphere() {};

StructuredCubedSphere::StructuredCubedSphere(const Field& tij, const Field& ghost) :
  tij_(tij), ghost_(ghost)
{

  Log::debug() << "StructuredCubedSphere bounds checking is set to "
                  + std::to_string(checkBounds) << std::endl;

  // Make array view.
  const auto tijView_ = array::make_view<idx_t, 2>(tij);

  // loop over tij and find min and max ij bounds.
  for (idx_t index = 0; index < tijView_.shape(0); ++index) {
    const size_t t = static_cast<size_t>(tijView_(index, Coordinates::T));
    const idx_t i = tijView_(index, Coordinates::I);
    const idx_t j = tijView_(index, Coordinates::J);

    ijBounds_[t].iBegin = std::min(i    , ijBounds_[t].iBegin);
    ijBounds_[t].jBegin = std::min(j    , ijBounds_[t].jBegin);
    ijBounds_[t].iEnd   = std::max(i + 1, ijBounds_[t].iEnd  );
    ijBounds_[t].jEnd   = std::max(j + 1, ijBounds_[t].jEnd  );

  }

  // Set tijToIdx vectors
  for (idx_t t = 0; t < 6; ++t) {

    // Set data array.
    const size_t vecSize = static_cast<size_t>((j_end(t) - j_begin(t))
                                           * (i_end(t) - i_begin(t)));
    tijToIdx_.push_back(std::vector<idx_t>(vecSize, invalid_index()));

  }

  // loop over ijt_ and set ijtToIdx
  for (idx_t index = 0; index < tijView_.shape(0); ++index) {
    const idx_t t = tijView_(index, Coordinates::T);
    const idx_t i = tijView_(index, Coordinates::I);
    const idx_t j = tijView_(index, Coordinates::J);

   tijToIdx_[static_cast<size_t>(t)][vecIndex(t, i, j)] = index;
  }
}

idx_t StructuredCubedSphere::i_begin(idx_t t) const {
  if (checkBounds) tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].iBegin;
}

idx_t StructuredCubedSphere::i_end(idx_t t) const {
  if (checkBounds) tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].iEnd;
}

idx_t StructuredCubedSphere::j_begin(idx_t t) const {
  if (checkBounds) tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].jBegin;
}

idx_t StructuredCubedSphere::j_end(idx_t t) const {
  if (checkBounds) tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].jEnd;
}

idx_t StructuredCubedSphere::index(idx_t t, idx_t i, idx_t j) const {

  if (checkBounds) iBoundsCheck(i, t);
  if (checkBounds) jBoundsCheck(j, t);

  return tijToIdx_[static_cast<size_t>(t)][vecIndex(t, i, j)];
}

Field StructuredCubedSphere::tij() const {
  return tij_;
}

void StructuredCubedSphere::tBoundsCheck(idx_t t) const {
  if (t < 0 || t > 5) throw_OutOfRange("t", t, 6);
}

void StructuredCubedSphere::jBoundsCheck(idx_t j, idx_t t) const {
  const idx_t jSize = j_end(t) - j_begin(t);
  j -= j_begin(t);
  if (j < 0 || j >= jSize) throw_OutOfRange("j - jMin", j, jSize);
}

void StructuredCubedSphere::iBoundsCheck(idx_t i, idx_t t) const {
  const idx_t iSize = i_end(t) - i_begin(t);
  i -= i_begin(t);
  if (i < 0 || i >= iSize) throw_OutOfRange("i - iMin", i, iSize);
}

size_t StructuredCubedSphere::vecIndex(idx_t t, idx_t i, idx_t j) const {
  return static_cast<size_t>((j - j_begin(t)) * (i_end(t) - i_begin(t)) + i - i_begin(t));
}



} // detail
} // functionspace
} // atlas
