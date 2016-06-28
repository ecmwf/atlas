/*
 * (C) Copyright 1996-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#ifndef atlas_grid_global_Global_h
#define atlas_grid_global_Global_h

#include "atlas/grid/Grid.h"


namespace atlas {
namespace grid {
namespace global {


/**
 * @brief Global Grid
 * This class is a base class for all grids that are global (on the sphere)
 */
class Global: public Grid {

  public:

    static std::string className();

    static std::string grid_type_str();

    Global( const Domain& dom=Domain::makeGlobal() );

};


}  // namespace global
}  // namespace grid
}  // namespace atlas


#endif

