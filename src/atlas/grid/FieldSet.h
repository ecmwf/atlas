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

#ifndef atlas_grid_FieldSet_H
#define atlas_grid_FieldSet_H

#include <vector>

#include "eckit/types/Types.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/memory/ScopedPtr.h"

#include "eckit/grib/GribHandle.h" ///< @todo this is to be removed

#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/grid/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace eckit
{
    class PathName;
    class DataHandle;
	namespace grib { class GribHandle; }
}

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/// Represents a set of fields
/// The order of the fields is kept

class FieldSet : public eckit::Owned {

public: // types

    typedef eckit::SharedPtr<FieldSet> Ptr;

public: // methods

	static FieldSet::Ptr create( );

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

    const Grid& grid() const { ASSERT( !empty() ); return *grid_; }
    Grid& grid() { ASSERT( !empty() ); return *grid_; }

    std::vector<std::string> field_names() const;

protected: // methods

	Field::Ptr create_field( eckit::grib::GribHandle& );

protected: // members

	Field::Vector fields_; ///< field handle storage

	Grid::Ptr     grid_;   ///< describes the grid (shared)

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
