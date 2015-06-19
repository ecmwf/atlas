/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_LocalGrid_H
#define atlas_LocalGrid_H

#include <cstddef>
#include <vector>

#include "eckit/memory/Builder.h"
#include "atlas/Grid.h"

namespace atlas {
namespace grids {

//------------------------------------------------------------------------------------------------------

class LocalGrid : public Grid {

  public: // methods

    LocalGrid(Grid *grid, const Domain& domain);

    virtual ~LocalGrid();

    virtual size_t npts() const;

    virtual void lonlat( std::vector<Point>& ) const;

    virtual std::string gridType() const;
    virtual GridSpec spec() const;

  private:  // methods

    virtual void print(std::ostream&) const;

    /// Human readable name
    /// May not be unique, especially when BoundBox is different
    virtual std::string shortName() const;

    /// Hash of the information this class unique
    virtual void hash(eckit::MD5&) const;

    /// Brute-force cropping of points by the passed domain
    /// To be used by the constructor only
    void cropPoints();

  private: // members

    std::auto_ptr<Grid> grid_;

    mutable std::string  shortName_;

    std::vector<Grid::Point> localPts_;

};

//------------------------------------------------------------------------------------------------------

} // namespace grids
} // namespace atlas

#endif
