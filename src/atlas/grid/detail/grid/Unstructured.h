/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Willem Deconinck
/// @author Tiago Quintino
/// @author Pedro Maciel
/// @date January 2015


#pragma once

#include <cstddef>
#include <vector>
#include "eckit/memory/ScopedPtr.h"
#include "atlas/grid/Grid.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {


class Unstructured : public Grid {

public:

  class Iterator: public Grid::Iterator {
  public:
    Iterator(const Unstructured& grid):
        grid_(grid),
        n_(0) {
    }

    virtual bool next(PointXY& xy) {
      if( n_ != grid_.points_->size() ) {
        xy = grid_.xy(n_++);
        return true;
      } else {
        return false;
      }
    }

  private:
    const Unstructured& grid_;
    size_t n_;
  };

public: // methods

    static std::string grid_type_str();

    static std::string className();

    /// Constructor taking a list of parameters
    Unstructured( const Config& );

    /// Constructor taking a list of points
    Unstructured( std::vector< Point >* pts );

    /// Constructor taking a mesh
    Unstructured( const mesh::Mesh& m );

    virtual ~Unstructured();

    virtual size_t npts() const;

    virtual void lonlat(std::vector< Point >&) const;

    virtual Spec spec() const;

    /// Human readable name
    virtual std::string shortName() const;
    virtual std::string gridType() const { return "unstructured"; }

    PointXY xy(size_t n) const { return PointXY( (*points_)[n] ); }

    virtual Iterator* iterator() const{ return new Iterator(*this); }


private: // methods

    virtual void print(std::ostream&) const;


    /// Hash of the lonlat array + BoundingBox
    virtual void hash(eckit::MD5&) const;

protected:

    /// Storage of coordinate points
    eckit::ScopedPtr< std::vector< Point > > points_;

    /// Cache for the shortName
    mutable std::string shortName_;

    /// Cache for the spec since may be quite heavy to compute
    mutable eckit::ScopedPtr<Grid::Spec> cached_spec_;

};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
