/*
 * (C) Copyright 1996-2013 ECMWF.
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

#ifndef atlas_grid_Field_H
#define atlas_grid_Field_H

#include <vector>

#include "eckit/types/Types.h"
#include "eckit/memory/NonCopyable.h"
#include "Grid.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

/// @todo this class will become polymorphic
/// @todo move it out of this header
class MetaData : public eckit::StringDict {
public: // methods

    typedef std::unique_ptr<MetaData> Ptr;

    MetaData();

};

//------------------------------------------------------------------------------------------------------

class FieldH : private eckit::NonCopyable {
public: // types

    typedef std::shared_ptr<FieldH> Ptr;
    typedef std::vector< FieldH::Ptr >  Vector;
    typedef std::vector< FieldH::Ptr >::const_iterator Iterator;
    typedef atlas::FieldT< double > Data;

public: // methods

    FieldH( Grid::Ptr, MetaData::Ptr, Data& );

    const Grid& grid() const { return *grid_; }
    Grid& grid() { return *grid_; }

    const MetaData& metadata() const { return *metadata_; }
    MetaData& metadata() { return *metadata_; }

    Data& data() { return data_; }
    const Data& data() const { return data_; }

protected:

    Grid::Ptr       grid_;      ///< describes the grid (shared)
    MetaData::Ptr   metadata_;  ///< describes the field (singly owned)
    Data&           data_;      ///< stores the field data (not owned)

};

//------------------------------------------------------------------------------------------------------

/// Represents a set of fields
/// The order of the fields is kept
class FieldSet : private eckit::NonCopyable {

public: // methods

    /// Takes ownership of the fields
    FieldSet( const FieldH::Vector& fields = FieldH::Vector() );

    const FieldH::Vector& fields() const { return fields_; }

    FieldH::Vector& fields() { return fields_; }

protected:

    FieldH::Vector fields_;

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
