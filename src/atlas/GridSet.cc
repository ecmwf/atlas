/*
 * (C) Copyright 1996-2015 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <eckit/exception/Exceptions.h>

#include "atlas/FieldSet.h"
#include "atlas/Grid.h"

#include "atlas/GridSet.h"


namespace atlas {


bool GridSet::has(const Grid& grid) const
{
  return has(grid.uid());
}


bool GridSet::has(const Grid::uid_t& uid) const
{
  for (std::vector< Grid::Ptr >::const_iterator g = grids_.begin(); g != grids_.end(); ++g )
  {
    if ( uid == (*g)->uid() )
      return true;
  }
  return false;
}


void GridSet::push_back(Grid::Ptr grid)
{
  if( !has(grid->uid()) )
    grids_.push_back(grid);
}


Grid::Ptr GridSet::grid(const Grid::uid_t& uid) const
{
  for( std::vector< Grid::Ptr >::const_iterator g = grids_.begin(); g != grids_.end(); ++g )
  {
    if( uid == (*g)->uid() )
      return (*g);
  }

  throw eckit::UserError("Gridset::grid: requested grid not found (uid=\""+uid+"\")");
}


}  // namespace atlas

