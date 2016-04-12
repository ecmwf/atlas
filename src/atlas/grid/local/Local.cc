/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <string>
#include "atlas/grid/local/Local.h"

namespace atlas {
namespace grid {
namespace local {

//------------------------------------------------------------------------------

std::string Local::className() { return "atlas.grid.local.Local"; }

//------------------------------------------------------------------------------

Local::Local(const Domain& d)
  : Grid(d)
{
}

//------------------------------------------------------------------------------

} // namespace local
} // namespace grid
} // namespace atlas
