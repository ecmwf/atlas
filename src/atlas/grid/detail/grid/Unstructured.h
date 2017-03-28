/*
 * (C) Copyright 1996-2017 ECMWF.
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
#include "atlas/grid/detail/grid/Grid.h"


namespace atlas {
namespace grid {
namespace detail {
namespace grid {


class Unstructured : public Grid {

public:

  class Iterator: public Grid::Iterator {
  public:
    Iterator(const Unstructured& grid, bool begin=true):
        grid_(grid),
        n_( begin ? 0 : grid_.points_->size() ) {
    }

    virtual bool next(PointXY& xy) {
      if( n_ != grid_.points_->size() ) {
        xy = grid_.xy(n_++);
        return true;
      } else {
        return false;
      }
    }

    virtual const PointXY operator *() const {
        return grid_.xy(n_);
    }

    virtual const Grid::Iterator& operator ++() {
        ++n_;
        return *this;
    }

    virtual bool operator ==(const Grid::Iterator &other) const {
      return n_ == static_cast<const Iterator&>(other).n_;
    }

    virtual bool operator !=(const Grid::Iterator &other) const {
      return n_ != static_cast<const Iterator&>(other).n_;
    }

  private:
    const Unstructured& grid_;
    size_t n_;
  };

public: // methods

    static std::string static_type() { return "unstructured"; }
    virtual std::string name() const;
    virtual std::string type() const { return static_type(); }

    /// Constructor taking a list of parameters
    Unstructured( const Config& );

    /// Constructor taking a list of points
    Unstructured( std::vector< PointXY >* pts );

    /// Constructor taking a mesh
    Unstructured( const mesh::Mesh& m );

    virtual ~Unstructured();

    virtual size_t npts() const;

    virtual Spec spec() const;


    PointXY xy(size_t n) const { return (*points_)[n]; }

    PointLonLat lonlat(size_t n) const { return projection_.lonlat((*points_)[n]); }

    virtual Iterator* begin() const{ return new Iterator(*this); }
    virtual Iterator* end()   const{ return new Iterator(*this,false); }

private: // methods

    virtual void print(std::ostream&) const;

    /// Hash of the lonlat array + BoundingBox
    virtual void hash(eckit::MD5&) const;

protected:

    /// Storage of coordinate points
    eckit::ScopedPtr< std::vector< PointXY > > points_;

    /// Cache for the shortName
    mutable std::string shortName_;

    /// Cache for the spec since may be quite heavy to compute
    mutable eckit::ScopedPtr<Grid::Spec> cached_spec_;

};


}  // namespace grid
}  // namespace detail
}  // namespace grid
}  // namespace atlas
