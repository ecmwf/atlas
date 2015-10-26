/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @date October 2015

#ifndef atlas_mesh_Elements_H
#define atlas_mesh_Elements_H

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas { class Array; }
namespace atlas { namespace mesh { class Nodes; } }
namespace atlas { namespace mesh { class ElementType; } }

namespace atlas {
namespace mesh {

/**
 * \brief Elements class that owns a collection of fields defined in Elements of the mesh
 */
class Elements : public eckit::Owned {

public: // methods

//-- Constructors

  Elements( const Nodes& nodes );
  virtual ~Elements();

//-- Accessors

// -- Modifiers

  void add( ElementType*, size_t nb_elements );
  void add( ElementType*, size_t nb_elements, const size_t connectivity[] );

// -- Data
private:
  const Nodes& nodes_;
  std::vector< eckit::SharedPtr<ElementType> > element_types_;
  std::vector< eckit::SharedPtr<Array> > connectivity_;
};

extern "C"
{
}

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
