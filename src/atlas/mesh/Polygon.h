/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Pedro Maciel
/// @author Willem Deconinck
/// @date September 2017

#pragma once

#include <vector>
#include "eckit/memory/Owned.h"
#include "atlas/grid/detail/partitioner/Polygon.h"
#include "atlas/library/config.h"
#include "atlas/mesh/detail/MeshImpl.h"

namespace atlas {
namespace mesh {

/**
 * \brief Polygon class that holds the boundary of a mesh partition
 */
class Polygon : public grid::detail::partitioner::Polygon, public eckit::Owned {

  typedef grid::detail::partitioner::Polygon Poly;

public: // methods

//-- Constructors

  /// @brief Construct "size" Polygon
  Polygon( const detail::MeshImpl& mesh, size_t halo );

//-- Accessors

  size_t halo() const { return halo_; }

  /// @brief Return the memory footprint of the Polygon
  size_t footprint() const;

  void outputPythonScript(const eckit::PathName&) const;

private:

  void print(std::ostream&) const;

  friend std::ostream& operator<<(std::ostream& s, const Polygon& p) {
    p.print(s);
    return s;
  }

private:

  const detail::MeshImpl& mesh_;
  size_t halo_;

};

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
