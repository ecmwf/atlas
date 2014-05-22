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
#include <memory>

#include "eckit/types/Types.h"
#include "eckit/memory/NonCopyable.h"

#include "atlas/mesh/Mesh.hpp"
#include "atlas/mesh/Field.hpp"
#include "atlas/mesh/Metadata.hpp"
#include "atlas/grid/Grid.h"

//------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//------------------------------------------------------------------------------------------------------

class FieldHandle : private eckit::NonCopyable {
public: // types

    typedef std::shared_ptr<FieldHandle> Ptr;
    typedef std::vector< FieldHandle::Ptr > Vector;

    typedef atlas::FieldT< double > Data;
    typedef atlas::Mesh Mesh;

public: // methods

    FieldHandle( Grid::Ptr, Data& );

    const Grid& grid() const { return *grid_; }
    Grid& grid() { return *grid_; }

    const Mesh& mesh() const { return *mesh_; }
    Mesh& mesh() { return *mesh_; }

    const Metadata& metadata() const { return data_.metadata(); }
    Metadata& metadata() { return data_.metadata(); }

    Data& data() { return data_; }
    const Data& data() const { return data_; }

protected: // members

    Grid::Ptr       grid_;      ///< describes the grid (shared)
    Mesh::Ptr       mesh_;      ///< mesh data structure (shared)
    Data&           data_;      ///< reference to the field data, not owned since it actually exists in the mesh_

};

//------------------------------------------------------------------------------------------------------

/// Represents a set of fields
/// The order of the fields is kept
class FieldSet : private eckit::NonCopyable {

public: // methods

    /// Takes ownership of the fields
    FieldSet( const FieldHandle::Vector& fields = FieldHandle::Vector() );

    const FieldHandle::Vector& fields() const { return fields_; }

    FieldHandle::Vector& fields() { return fields_; }

protected:

    FieldHandle::Vector fields_;

};

//------------------------------------------------------------------------------------------------------

} // namespace grid
} // namespace atlas

#endif
