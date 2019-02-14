/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

/// @author Pedro Maciel
/// @author Willem Deconinck
/// @date September 2017

#pragma once

#include <vector>

#include "atlas/util/Object.h"

#include "atlas/library/config.h"
#include "atlas/util/Config.h"
#include "atlas/util/Polygon.h"

namespace atlas {
namespace mesh {
namespace detail {
class MeshImpl;
}
}  // namespace mesh
}  // namespace atlas

namespace atlas {
namespace mesh {

/**
 * @brief Polygon class that holds the boundary of a mesh partition
 */
class PartitionPolygon : public util::Polygon, public util::Object {
public:  // methods
    //-- Constructors

    /// @brief Construct "size" polygon
    PartitionPolygon( const detail::MeshImpl& mesh, idx_t halo );

    //-- Accessors

    idx_t halo() const { return halo_; }

    /// @brief Return the memory footprint of the Polygon
    size_t footprint() const;

    void outputPythonScript( const eckit::PathName&, const eckit::Configuration& = util::NoConfig() ) const;

private:
    void print( std::ostream& ) const;

    friend std::ostream& operator<<( std::ostream& s, const PartitionPolygon& p ) {
        p.print( s );
        return s;
    }

private:
    const detail::MeshImpl& mesh_;
    idx_t halo_;
};

//------------------------------------------------------------------------------------------------------

}  // namespace mesh
}  // namespace atlas
