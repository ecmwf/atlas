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


#include "atlas/Grid.h"


namespace atlas {


/**
 * @brief GridSet contains the Grid objects whithin a consistent context (p.e. the execution of MIR)
 */
class GridSet {

public: // methods

  bool has( const Grid& ) const;
  bool has( const Grid::uid_t& ) const;

  size_t size() const { return grids_.size(); }

  void clear()              { grids_.clear(); }

  const Grid& operator[](size_t i) const { return *(grids_.at(i)); }
  Grid&       operator[](size_t i)       { return *(grids_.at(i)); }

  /**
   * @brief insert a grid in the container, if it is unique to it
   * @param grid grid to insert
   */
  void push_back(Grid::Ptr grid);

  /**
   * @return grid based on unique identifier (a string)
   * @param uid grid identifier
   */
  Grid::Ptr grid(const Grid::uid_t& uid) const;

private: // data

  /// Container for Grid's
  std::vector< Grid::Ptr > grids_;

};


} // namespace atlas


#endif
