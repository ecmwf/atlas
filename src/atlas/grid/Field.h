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

#ifndef eckit_grid_Field_H
#define eckit_grid_Field_H

#include <vector>

#include "eckit/types/Types.h"
#include "eckit/memory/NonCopyable.h"
#include "Grid.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {


//-----------------------------------------------------------------------------

/// Represents a Field which store one type of data parameter
class Field : private eckit::NonCopyable {
public: // types

    class MetaData : public eckit::StringDict {
    public: // methods

        MetaData();

    };

    typedef std::vector< Field* >  Vector;
    typedef std::vector< Field* >::const_iterator Iterator;
    typedef std::vector< double >  Data;

public: // methods

    Field( Grid* grid, MetaData* metadata, std::vector<double>* data );

    ~Field();

    const MetaData& metadata() const { return *metadata_; }

    Data& data() { return *data_; }
    const Data& data() const { return *data_; }

    /// @returns the grid reference
    const Grid& grid() const { return *grid_; }


protected:

    /// @todo make the grid a shared pointer

    Grid*       grid_;      ///< describes the grid

    MetaData*   metadata_;  ///< describes the field

    Data*       data_;      ///< stores the field data

};

//-----------------------------------------------------------------------------

/// Represents a set of fields
/// The order of the fields is kept
class FieldSet : private eckit::NonCopyable {

public: // methods

    /// Takes ownership of the fields
    FieldSet( const Field::Vector& fields = Field::Vector() );

    ~FieldSet();
    
    const Field::Vector& fields() const { return fields_; }

    Field::Vector& fields() { return fields_; }

protected:

    Field::Vector fields_;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
