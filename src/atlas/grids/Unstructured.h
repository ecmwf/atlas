/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_Unstructured_H
#define atlas_Unstructured_H

/// @author Tiago Quintino
/// @author Pedro Maciel
/// @date January 2015

#include <cstddef>
#include <vector>

#include "eckit/memory/ScopedPtr.h"
#include "atlas/Grid.h"


namespace atlas {
namespace grids {


class Unstructured : public Grid {

public: // methods

  static std::string grid_type_str() { return "unstructured"; }

  static std::string className() { return "atlas.grid.Unstructured"; }

  /// Constructor taking a list of parameters
  Unstructured(const eckit::Params& p);

  /// Constructor taking a list of points
  Unstructured(std::vector< Point >* pts);

  virtual ~Unstructured();

  virtual BoundBox bounding_box() const;

  virtual size_t npts() const;

  virtual void lonlat(double[]) const;
  virtual void lonlat(std::vector< Point >&) const;
  virtual void lonlat(std::vector< double >&) const;

  virtual std::string grid_type() const { return grid_type_str(); }

  virtual GridSpec spec() const;

private: // methods

  /// Human readable name
  virtual std::string shortName() const;

  /// Hash of the latlon array + BoundBox
  virtual void hash(eckit::MD5&) const;

protected:

  eckit::ScopedPtr< std::vector< Point > > points_;  ///< storage of coordinate points

  BoundBox bound_box_;                               ///< bounding box for the domain

  mutable std::string shortName_;      ///< cache for the shortName

  mutable eckit::ScopedPtr<GridSpec> cachedGridSpec_;  ///< cache for the GridSpec since may be quite heavy to compute

};


} // namespace grids
} // namespace atlas
#endif
