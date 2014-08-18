/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "atlas/grid/GridSpecParams.h"

using namespace eckit;
using namespace std;

namespace atlas {
namespace grid {


GridSpecParams::GridSpecParams( const GridSpec& p) : spec_(p)
{
   ASSERT(!p.empty());
   ASSERT(!spec_.empty());
}

Params::value_t GridSpecParams::get( const key_t& key ) const
{
#ifdef DEBUG_ME
   Params::value_t v = spec_.get(key);
   if (v.isNil()) {
       std::cout << "GridSpecParams::get failed for key " << key << "  Opps" << std::endl;
   }
   return v;
#endif

   return spec_.get(key);
}

void GridSpecParams::print(std::ostream& s) const
{
   s << spec_;
}

// ------------------------------------------------------------------

} // namespace grid
} // namespace atlas

