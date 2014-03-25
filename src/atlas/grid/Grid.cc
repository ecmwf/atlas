/*
 * (C) Copyright 1996-2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/log/Log.h"

#include "eckit/grid/Grid.h"

//-----------------------------------------------------------------------------

namespace eckit {
namespace grid {

//-----------------------------------------------------------------------------

Grid::Grid()
{
    Log::info() << "Build a Grid" << std::endl;
}

Grid::~Grid()
{
    Log::info() << "Destroy a Grid" << std::endl;
}

std::string Grid::hash() const
{
    NOTIMP;
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit
