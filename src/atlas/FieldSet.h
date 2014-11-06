/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Peter Bispham
/// @author Tiago Quintino
/// @date Oct 2013

#ifndef atlas_FieldSet_H
#define atlas_FieldSet_H

#include <vector>

#include <eckit/types/Types.h>
#include <eckit/memory/Owned.h>
#include <eckit/memory/SharedPtr.h>
#include <eckit/memory/ScopedPtr.h>

#include "atlas/atlas_config.h"

#ifdef ECKIT_HAVE_GRIB
  #include <eckit/grib/GribHandle.h> ///< @todo this is to be removed
#endif

#include "atlas/Mesh.h"
#include "atlas/Field.h"
#include "atlas/Metadata.h"
#include "atlas/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace eckit
{
  class PathName;
  class DataHandle;

  namespace grib { class GribHandle; }
}

namespace atlas {


//------------------------------------------------------------------------------------------------------

/// Represents a set of fields
/// The order of the fields is kept

class FieldSet : public eckit::Owned {

public: // types

    typedef eckit::SharedPtr<FieldSet> Ptr;

public: // methods

	/// Constructs a field set from a file (e.g. a GRIB file )
    FieldSet( const eckit::PathName& );

    /// @todo Constructor for a FieldSet from a buffer
    FieldSet( const eckit::Buffer& );

    /// @todo Constructor for a FieldSet from a DataHandle
    //  FieldSet( const eckit::DataHandle& );

    /// Constructs a field set with n fields from a Grid
    FieldSet( const Grid::Ptr grid, std::vector<std::string> nfields );

    /// Constructs a FielSet from predefined fields
    /// Takes ownership of the fields
	FieldSet( const Field::Vector& fields );

	const Field& operator[]( const size_t& i ) const { ASSERT(i<size()); return *fields_[i]; }
	Field& operator[]( const size_t& i )             { ASSERT(i<size()); return *fields_[i]; }

	const Field::Vector& fields() const { return fields_; }
	Field::Vector& fields() { return fields_; }

    size_t size() const { return fields_.size(); }
    bool empty() const { return ! fields_.size(); }

	const Grid& grid() const { ASSERT( !empty() ); return fields_[0]->grid(); }
	Grid& grid() { ASSERT( !empty() ); return fields_[0]->grid(); }

    std::vector<std::string> field_names() const;

private: // methods

  Field::Ptr create_field( eckit::grib::GribHandle& );

	bool checkConsistency() const;

protected: // members

	Field::Vector fields_; ///< field handle storage

	Grid::Ptr grid_;

};

//------------------------------------------------------------------------------------------------------


} // namespace atlas

#endif
