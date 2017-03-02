/*
 * (C) Copyright 1996-2017 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#pragma once

#include "atlas/grid/detail/grid/Grid.h"

//---------------------------------------------------------------------------------------------------------------------

namespace atlas {
namespace grid {

//---------------------------------------------------------------------------------------------------------------------

class Iterator {

public:

    Iterator( detail::grid::Grid::Iterator* iterator ):
        iterator_(iterator) {
    }

    bool next( PointXY& xy ) { return iterator_->next(xy); }

private:

    std::unique_ptr<detail::grid::Grid::Iterator> iterator_;
};

//---------------------------------------------------------------------------------------------------------------------

}
}
