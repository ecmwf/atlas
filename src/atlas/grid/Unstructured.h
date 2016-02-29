/*
 * (C) Copyright 1996-2016 ECMWF.
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
#include "atlas/grid/Grid.h"


namespace atlas {
namespace grid {


class Unstructured : public Grid {

public: // methods

  static std::string grid_type_str() { return "unstructured"; }

  static std::string className() { return "atlas.grid.Unstructured"; }

  /// Constructor taking a list of parameters
  Unstructured(const eckit::Parametrisation& p);

  /// Constructor taking a list of points
  Unstructured(std::vector< Point >* pts);

  /// Constructor taking a mesh
  Unstructured(const mesh::Mesh& m);

  virtual ~Unstructured();

  virtual BoundBox boundingBox() const;

  virtual size_t npts() const;

  virtual void lonlat(std::vector< Point >&) const;

  virtual std::string gridType() const { return grid_type_str(); }

  virtual eckit::Properties spec() const;

private: // methods

  virtual void print(std::ostream&) const;

  /// Human readable name
  virtual std::string shortName() const;

  /// Hash of the lonlat array + BoundBox
  virtual void hash(eckit::MD5&) const;

protected:

  eckit::ScopedPtr< std::vector< Point > > points_;  ///< storage of coordinate points

  BoundBox bound_box_;                               ///< bounding box for the domain

  mutable std::string shortName_;      ///< cache for the shortName

  mutable eckit::ScopedPtr<eckit::Properties> cached_spec_;  ///< cache for the spec since may be quite heavy to compute

};


} // namespace grid
} // namespace atlas
#endif
