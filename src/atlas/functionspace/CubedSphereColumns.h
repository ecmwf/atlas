/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/CellColumns.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace/detail/CubedSphereColumnsImpl.h"
#include "atlas/meshgenerator/detail/cubedsphere/CubedSphereUtility.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
class Mesh;
class CellColumns;
class NodeColumns;
}

namespace atlas {
namespace functionspace {

/// Extend NodeColumns and CellColumns so that they can exploit CubedSphere structure.
template <typename BaseFunctionSpace>
class CubedSphereColumns : public BaseFunctionSpace {

public:

  /// Constructors.
  CubedSphereColumns();
  CubedSphereColumns(const FunctionSpace& functionSpace);
  CubedSphereColumns(const Mesh& mesh, const eckit::Configuration& configuration);
  CubedSphereColumns(const Mesh& mesh);

  /// Invalid index.
  static constexpr idx_t invalid_index() {return -1;}

  /// i lower bound for tile t (including halo)
  idx_t i_begin(idx_t t) const;
  /// i lower bound for tile t (including halo)
  idx_t i_end(idx_t t) const;

  /// j lower bound for tile t (including halo)
  idx_t j_begin(idx_t t) const;
  /// j lower bound for tile t (including halo)
  idx_t j_end(idx_t t) const;

  /// Return array_view index for (t, i, j).
  idx_t index(idx_t t, idx_t i, idx_t j) const;

  /// Return tij field.
  Field tij() const;

  /// Apply functor f(index, t, i, j) to all valid (t, i, j) triplets.
  /// (use "#pragma omp atomic" to update array elements other than [index].)
  template<typename functor>
  void for_each(const functor& f, bool include_halo = false)  const {

    using namespace meshgenerator::detail::cubedsphere;

    // make array views.
    const auto tijView_ =
      array::make_view<idx_t, 2>( cubedSphereColumnsHandle_.get()->tij() );
    const auto ghostView_ =
      array::make_view<idx_t, 1>( cubedSphereColumnsHandle_.get()->ghost() );

    // Loop over elements.
    for (idx_t index = 0; index < tijView_.shape(0); ++index) {

      if (!include_halo and ghostView_(index)) continue;

      const idx_t t = tijView_(index, Coordinates::T);
      const idx_t i = tijView_(index, Coordinates::I);
      const idx_t j = tijView_(index, Coordinates::J);

      f(index, t, i, j);

    }
  }

private:

  // Object hanldle for CubedSphereColumnsImpl
  util::ObjectHandle<detail::CubedSphereColumnsImpl> cubedSphereColumnsHandle_;

};

using CubedSphereNodeColumns = CubedSphereColumns<NodeColumns>;
using CubedSphereCellColumns = CubedSphereColumns<CellColumns>;

} // namespace functionspace
} // namespace atlas
