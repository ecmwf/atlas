/*
 * (C) Copyright 1996-2014 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef atlas_grid_GridSpecParams_H
#define atlas_grid_GridSpecParams_H

#include "eckit/value/Params.h"
#include "atlas/grid/GridSpec.h"

namespace atlas {
namespace grid {

class GridSpecParams : public eckit::Params {

public: // methods

   GridSpecParams( const GridSpec& );

   virtual value_t get( const key_t& key ) const;

protected: // methods

   virtual void print(std::ostream& s) const;

private: // members

   GridSpec spec_;
};

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace atlas


#endif
