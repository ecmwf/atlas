#ifndef atlas_grid_Unstructured_H
#define atlas_grid_Unstructured_H
/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Tiago Quintino
/// @date April 2013


#include <cstddef>
#include <vector>

#include "eckit/memory/ScopedPtr.h"

#include "atlas/grid/Grid.h"
#include "atlas/grid/GridSpec.h"

//-----------------------------------------------------------------------------

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

class Unstructured : public Grid {

public: // methods

    /// @warning temporary constructor taking a list of points
    Unstructured( std::vector< Point >* pts, const std::string& hash );

    virtual ~Unstructured();

    virtual std::string hash() const;

    virtual BoundBox boundingBox() const;

    virtual size_t nPoints() const;

    virtual void coordinates( Grid::Coords & ) const;

    virtual const GridSpec& spec() const { return the_grid_spec_ ;}

    /// @deprecated will be removed soon as it exposes the inner storage of the coordinates
    virtual const std::vector<Point>& coordinates() const { return *points_; }

protected:

    eckit::ScopedPtr< std::vector< Point > > points_; ///< storage of coordinate points

    BoundBox bound_box_;              ///< bounding box for the domain

    std::string hash_;

    GridSpec   the_grid_spec_;       ///< unique description of Grid
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

#endif
