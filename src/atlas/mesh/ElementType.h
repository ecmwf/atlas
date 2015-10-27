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

#ifndef atlas_mesh_ElementType_H
#define atlas_mesh_ElementType_H

#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"

namespace atlas {
namespace mesh {

/**
 * \brief ElementType class (abstract) that provides access to geometric information of an element
 */
class ElementType : public eckit::Owned {

public: // methods

  static ElementType* create( const std::string& );

//-- Constructors

  ElementType();
  ~ElementType() = 0;

//-- Accessors

  virtual const std::string& name() const = 0;
  virtual size_t nb_nodes() const = 0;
  virtual size_t nb_edges() const = 0;

};

extern "C"
{
}

//------------------------------------------------------------------------------------------------------

} // namespace mesh
} // namespace atlas

#endif
