/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/grid/Iterator.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

IterateXY::iterator IterateXY::begin() const {
    return grid_.xy_begin();
}

IterateXY::iterator IterateXY::end() const {
    return grid_.xy_end();
}

IterateLonLat::iterator IterateLonLat::begin() const {
    return grid_.lonlat_begin();
}

IterateLonLat::iterator IterateLonLat::end() const {
    return grid_.lonlat_end();
}

}  // namespace grid
}  // namespace atlas
