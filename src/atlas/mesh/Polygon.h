/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date August 2015

#pragma once

#include <vector>
#include "eckit/memory/Owned.h"
#include "atlas/mesh/detail/MeshImpl.h"
#include "atlas/library/config.h"

namespace atlas {
namespace mesh {

/**
 * \brief Polygon class that holds the boundary of a mesh partition
 */
class Polygon : public eckit::Owned {
public:

  using const_iterator = std::vector<idx_t>::const_iterator;

public: // methods

//-- Constructors

  /// @brief Construct "size" Polygon
  Polygon( const detail::MeshImpl& mesh, size_t halo );

//-- Accessors
  size_t size() const { return node_list_.size(); }

  size_t halo() const { return halo_; }

  /// @brief Return the memory footprint of the Polygon
  size_t footprint() const;

  void outputPythonScript(const eckit::PathName&) const;

  const_iterator begin() const { return node_list_.begin(); }
  const_iterator end()   const { return node_list_.end();   }

private:

  void print(std::ostream&) const;

  friend std::ostream& operator<<(std::ostream& s, const Polygon& p) {
    p.print(s);
    return s;
  }

private:

  const detail::MeshImpl& mesh_;
  size_t halo_;
  std::vector<idx_t> node_list_;

};

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas
