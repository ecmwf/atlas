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
#include "atlas/grid/global/Global.h"


namespace atlas {
namespace grid {
namespace global {


std::string Global::className() {
    return "atlas.grid.global.Global";
}


std::string Global::grid_type_str() {
    return "global";
}


Global::Global(const Domain& dom) :
    Grid(dom) {
}


}  // namespace global
}  // namespace grid
}  // namespace atlas

