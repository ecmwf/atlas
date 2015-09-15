/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

/// @author Peter Bispham
/// @author Tiago Quintino
/// @author Pedro Maciel
/// @date Jan 2015

#ifndef atlas_Grid_H
#define atlas_Grid_H

#include <cstddef>
#include <vector>
#include <cmath>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "eckit/value/Properties.h"
#include "eckit/geometry/Point2.h"
#include "eckit/memory/Builder.h"
#include "eckit/memory/Owned.h"
#include "eckit/memory/SharedPtr.h"
#include "eckit/utils/MD5.h"

#include "atlas/BoundBox.h"
#include "atlas/util/ObjectRegistry.h"
#include "atlas/Domain.h"
#include "atlas/Config.h"

namespace atlas {

class Mesh;

//------------------------------------------------------------------------------------------------------

class Grid : public eckit::Owned, public util::Registered<Grid> {

 public:  // types

  typedef eckit::BuilderT1<Grid> builder_t;
  typedef const eckit::Parametrisation& ARG1;

  typedef eckit::SharedPtr<Grid> Ptr;

  typedef eckit::geometry::LLPoint2 Point; // must be sizeof(double)*2
  typedef std::string uid_t;

 public:  // methods

  static std::string className() { return "atlas.Grid"; }

  static Grid* create(const eckit::Parametrisation&);
  static Grid* create(const Grid::uid_t& shortName);

  /// Default constructor builds a Grid that is Global
  Grid();

  /// Constructor takes a domain for this Grid
  Grid( const Domain& );

  virtual ~Grid();

  Ptr self() { return Ptr(this); } ///< @todo not necessary?

  /// Human readable name (may not be unique)
  virtual std::string shortName() const = 0;

  /// Unique grid id
  /// Computed from the shortName and the hash
  uid_t uniqueId() const;

  /// Adds to the MD5 the information that makes this Grid unique
  virtual void hash(eckit::MD5&) const = 0;

  /// @returns the hash of the information that makes this Grid unique
  eckit::MD5::digest_t hash() const;

  ///< @todo not necessary?

  /// @return bounding box of all lonlat points,
  /// @note this function will compute brute-force the minmax lonlat values over all the points
  virtual BoundBox boundingBox() const;

  /// @return area which contains the grid
  virtual const Domain& domain() const;

  /// @return number of grid points
  /// @note This methods should have constant access time, if necessary derived
  //        classes should compute it at construction
  virtual size_t npts() const = 0;

  /// Fill provided parameter with grid points, as (lon,lat) values
  /// @post resizes the vector
  virtual void lonlat(std::vector<Point>&) const = 0;

  /// Fills the provided vector with the (lon,lat) values
  /// @post resizes the vector
  void fillLonLat(std::vector<double>&) const;

  /// Fills the provided array with the (lon,lat) values
  /// @note Assumes that the input array has been allocated with correct size
  /// @param array is an array already allocated with enough size to store all the latlon values
  /// @param arraySize is the size of the array
  void fillLonLat(double array[], size_t arraySize) const;

  virtual std::string gridType() const = 0;

  virtual eckit::Properties spec() const = 0;

  virtual bool same(const Grid&) const;

  /// @TODO: eventually remove the Mesh from the Grid

  void set_mesh(const Mesh& mesh);
  Mesh& mesh() const;

protected:  // methods

  /// Fill provided memory buffer with the grid points, as (lon,lat) values
  /// This implementation in the base Grid class is not optimal as it incurs in double copy
  /// Derived classes should reimplement more optimised versions.
  ///
  /// @note Assumes that the input buffer has been allocated with correct size,
  ///       possibly from calling method npts()
  ///
  /// @param array to be filled in with the (lon,lat) values
  /// @param size number of doubles in array
  ///
  /// @return the size of bytes copyied in
  virtual size_t copyLonLatMemory(double* pts, size_t size) const;

  virtual void print(std::ostream&) const = 0;

private:  // methods

  friend std::ostream& operator<<(std::ostream& s, const Grid& p) {
      p.print(s);
      return s;
  }

protected:  // members

  atlas::Domain domain_;

private:  // members

  mutable eckit::SharedPtr<Mesh> mesh_; ///< @todo to be removed

  mutable uid_t                  uid_;  ///< cache the unique ID
  mutable eckit::MD5::digest_t   hash_; ///< cache the hash

};

//------------------------------------------------------------------------------------------------------

}  // namespace atlas

#endif
