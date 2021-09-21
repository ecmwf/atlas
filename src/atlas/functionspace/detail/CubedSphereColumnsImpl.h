/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/util/Object.h"

/// Extra cubed sphere functionality for NodeColumns and CellColumns function
/// spaces.

namespace atlas {
class Field;
}

namespace atlas {
namespace functionspace {
namespace detail {

using namespace meshgenerator::detail::cubedsphere;

class CubedSphereColumnsImpl : public util::Object
{

public:

  CubedSphereColumnsImpl() = default;
  CubedSphereColumnsImpl(const Field& tij, const Field& ghost);

  /// Invalid index.
  static constexpr idx_t invalid_index() {return -1;}

  /// i lower bound for tile t (including halo)
  idx_t i_begin(idx_t) const;
  /// i lower bound for tile t (including halo)
  idx_t i_end(idx_t t) const;

  /// j lower bound for tile t (including halo)
  idx_t j_begin(idx_t t) const;
  /// j lower bound for tile t (including halo)
  idx_t j_end(idx_t t) const;

  /// Return array_view index for (i, j, t).
  idx_t index(idx_t t, idx_t i, idx_t j) const;

  /// Return ijt field.
  Field tij() const;

  /// Return ghost/halo field.
  Field ghost() const;

private:

  // Bounding box struct.
  struct BoundingBox {
      idx_t iBegin{std::numeric_limits<idx_t>::max()};
      idx_t iEnd{std::numeric_limits<idx_t>::min()};
      idx_t jBegin{std::numeric_limits<idx_t>::max()};
      idx_t jEnd{std::numeric_limits<idx_t>::min()};
  };

  // Bounds checking.
  void tBoundsCheck(idx_t t) const;
  void iBoundsCheck(idx_t i, idx_t t) const;
  void jBoundsCheck(idx_t j, idx_t t) const;

  // Row-major index of ijtToidx vectors.
  size_t vecIndex(idx_t i, idx_t j, idx_t t) const;

  // Index storage vectors.
  std::vector<std::vector<idx_t>> tijToIdx_;

  // Array of bounding boxes.
  std::array<BoundingBox, 6> ijBounds_;

  // ijt field.
  Field tij_;

  // ghost field.
  Field ghost_;

};

} // detail
} // functionspace
} // atlas
