/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date October 2015

#pragma once

#include "atlas/mesh/elementtypes/Classification.h"

namespace atlas {
namespace mesh {
namespace elementtypes {

//------------------------------------------------------------------------------------------------------

class Line : public Edge
{
public:
  enum { VERTICES = 2 };
  enum { FACETS   = VERTICES };
  enum { RIDGES   = 0 };
  virtual ~Line() {}
  virtual bool parametric() const { return true; }
  virtual bool simplex() const { return false; }
  virtual size_t nb_vertices() const { return VERTICES; }
  virtual size_t nb_edges()    const { return EDGES;    }
  virtual size_t nb_nodes()    const { return VERTICES; }
  virtual size_t nb_facets()   const { return FACETS;   }
  virtual const std::string& name() const { static std::string s("Line"); return s; }
};

//------------------------------------------------------------------------------------------------------

} // namespace elementtypes

namespace temporary {
  class Line : public elementtypes::Line {
    public:
    [[deprecated("Use 'atlas::mesh::ElementType::create(\"Line\")' instead")]] Line() {}
  };
}

} // namespace mesh
} // namespace atlas
