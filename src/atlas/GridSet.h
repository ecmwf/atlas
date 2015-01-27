/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_GridSet_H
#define atlas_GridSet_H


#include "eckit/memory/SharedPtr.h"


// forward declarations
namespace atlas {
  class FieldSet;
  class Grid;
}

namespace atlas {


/**
 * @brief GridSet contains the Grid objects whithin a consistent context (p.e. the execution of MIR)
 */
class GridSet /*: private std::vector< eckit::SharedPtr< Grid > >*/ {
#if 0
public:
  // vector-like functionality (expose some internals)

  using vector::iterator;
  using vector::const_iterator;
  using vector::reverse_iterator;
  using vector::const_reverse_iterator;
  using vector::begin;
  using vector::end;

  using vector::size;
  using vector::clear;
  using vector::pop_back;
  using vector::swap;

  using vector::operator[];

  using vector::push_back;
  void push_back(Grid* grid) { push_back(eckit::SharedPtr< Grid >(grid)); }

#else
public:
  // vector-like functionality (mimic syntax)

  size_t size() const { return grids_.size(); }
  void clear()              { grids_.clear(); }
  void pop_back()           { grids_.pop_back(); }
  void swap(GridSet& other) { grids_.swap(other.grids_); }

  const Grid& operator[](size_t i) const { return *(grids_.at(i)); }
        Grid& operator[](size_t i)       { return *(grids_.at(i)); }

  void push_back(const eckit::SharedPtr< Grid > grid) {grids_.push_back(grid); }
#endif

public:
  // specific functionality

  /**
   * @return grid based on unique identifier (a string)
   * @param uid grid identifier
   */
  const Grid& grid(const std::string& uid) const;

  /**
   * @return grid based on unique identifier (a string)
   * @param uid grid identifier
   */
  Grid& grid(const std::string& uid);

private:
  // data

  /// Container for Grid's
  std::vector< eckit::SharedPtr< Grid > > grids_;

};


} // namespace atlas

#endif

