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


const Grid& GridSet::grid(const std::string &uid) const
{
#if 0
  for (const_iterator g=begin(); g!=end(); ++g)
    if (uid==(*g)->uid())
      return (*(*g));
  throw eckit::UserError("Gridset::grid: requested grid not found (uid=\""+uid+"\")");
  return (*(*begin()));
#else
  for (std::vector< eckit::SharedPtr< Grid > >::const_iterator g=grids_.begin(); g!=grids_.end(); ++g)
    if (uid==(*g)->uid())
      return (*(*g));
  throw eckit::UserError("Gridset::grid: requested grid not found (uid=\""+uid+"\")");
  (*(*grids_.begin()));
#endif
}


Grid& GridSet::grid(const std::string &uid)
{
  return const_cast< Grid& >(grid(uid));
}


}  // namespace atlas

