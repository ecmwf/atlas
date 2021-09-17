/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "CubedSphereUtility.h"

namespace atlas {
namespace functionspace {
namespace detail {

CubedSphereUtility::~CubedSphereUtility() {};

CubedSphereUtility::CubedSphereUtility(const Field& ijt, const Field& ghost) :
  ijt_(ijt), ghost_(ghost)
{

  // Make array view.
  const auto ijtView_ = array::make_view<idx_t, 2>(ijt);

  // loop over ijt_ and find min and max ij bounds.
  for (idx_t index = 0; index < ijtView_.shape(0); ++index) {
    const auto i = ijtView_(index, Coordinates::I);
    const auto j = ijtView_(index, Coordinates::J);
    const auto t = static_cast<size_t>(ijtView_(index, Coordinates::T));

    ijBounds_[t].iBegin = std::min(i    , ijBounds_[t].iBegin);
    ijBounds_[t].jBegin = std::min(j    , ijBounds_[t].jBegin);
    ijBounds_[t].iEnd   = std::max(i + 1, ijBounds_[t].iEnd  );
    ijBounds_[t].jEnd   = std::max(j + 1, ijBounds_[t].jEnd  );

  }

  // Set ijtToIdx vectors
  for (idx_t t = 0; t < 6; ++t) {

    // Set data array.
    const auto vecSize = (j_end(t) - j_begin(t))
                       * (i_end(t) - i_begin(t));
    ijtToIdx_.push_back(std::vector<idx_t>(static_cast<size_t>(vecSize), invalid_index()));

  }

  // loop over ijt_ and set ijtToIdx
  for (idx_t index = 0; index < ijtView_.shape(0); ++index) {
    const auto i = ijtView_(index, Coordinates::I);
    const auto j = ijtView_(index, Coordinates::J);
    const auto t = ijtView_(index, Coordinates::T);

   ijtToIdx_[static_cast<size_t>(t)][vecIndex(i, j, t)] = index;
  }
}

idx_t CubedSphereUtility::i_begin(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].iBegin;
}

idx_t CubedSphereUtility::i_end(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].iEnd;
}

idx_t CubedSphereUtility::j_begin(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].jBegin;
}

idx_t CubedSphereUtility::j_end(idx_t t) const {
  tBoundsCheck(t);
  return ijBounds_[static_cast<size_t>(t)].jEnd;
}

idx_t CubedSphereUtility::index(idx_t i, idx_t j, idx_t t) const {

  // Check bounds.
  iBoundsCheck(i, t);
  jBoundsCheck(j, t);

  return ijtToIdx_[static_cast<size_t>(t)][vecIndex(i, j, t)];
}

Field CubedSphereUtility::ijt() const {
  return ijt_;
}

void CubedSphereUtility::tBoundsCheck(idx_t t) const {
  if (t < 0 or t > 5) throw_OutOfRange("t", t, 6);
}

void CubedSphereUtility::jBoundsCheck(idx_t j, idx_t t) const {
  const auto jSize = j_end(t) - j_begin(t);
  j -= j_begin(t);
  if (j < 0 or j >= jSize) throw_OutOfRange("j - jMin", j, jSize);
}

void CubedSphereUtility::iBoundsCheck(idx_t i, idx_t t) const {
  const auto iSize = i_end(t) - i_begin(t);
  i -= i_begin(t);
  if (i < 0 or i >= iSize) throw_OutOfRange("i - iMin", i, iSize);
}

size_t CubedSphereUtility::vecIndex(idx_t i, idx_t j, idx_t t) const {
  return static_cast<size_t>((j - j_begin(t)) * (i_end(t) - i_begin(t))
               + i - i_begin(t));
}



} // detail
} // functionspace
} // atlas
