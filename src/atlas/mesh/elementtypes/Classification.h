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

#include "atlas/mesh/ElementType.h"

namespace atlas {
namespace mesh {
namespace elementtypes {

//-------------------------------------------------------------------

class Volume : public ElementType
{
public:
  enum { DIMENSIONALITY=3 };
  virtual size_t dimensionality() const { return DIMENSIONALITY; }
};

//-------------------------------------------------------------------

class Face : public ElementType
{
public:
  enum { DIMENSIONALITY=2 };
  enum { FACES=1 };
  virtual size_t dimensionality() const { return DIMENSIONALITY; }
  virtual size_t nb_faces() const { return FACES; }
};

//-------------------------------------------------------------------

class Edge: public ElementType
{
public:
  enum { DIMENSIONALITY=1 };
  enum { FACES=0 };
  enum { EDGES=1 };
  virtual size_t dimensionality() const { return DIMENSIONALITY; }
  virtual size_t nb_faces()   const { return FACES; }
  virtual size_t nb_edges()   const { return EDGES; }
};

//-------------------------------------------------------------------

class Vertex: public ElementType
{
public:
  enum { DIMENSIONALITY=0 };
  enum { FACES=0 };
  enum { EDGES=0 };
  enum { VERTICES=1 };
  virtual size_t dimensionality() const { return DIMENSIONALITY; }
  virtual size_t nb_faces()    const { return FACES;    }
  virtual size_t nb_edges()    const { return EDGES;    }
  virtual size_t nb_vertices() const { return VERTICES; }
};

//-------------------------------------------------------------------

} // namespace elementtypes
} // namespace mesh
} // namespace atlas
