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

#include "atlas/grid/Field.h"

using namespace eckit;

namespace atlas {
namespace grid {

//-----------------------------------------------------------------------------

MetaData::MetaData()
{
}

//-----------------------------------------------------------------------------

FieldH::FieldH( Grid::Ptr g, MetaData::Ptr md, Data& d ) :
    grid_(g),
    metadata_(std::move(md)),
    data_(d)
{
    ASSERT(grid_);
    ASSERT(metadata_);
}

//-----------------------------------------------------------------------------

FieldSet::FieldSet(const FieldH::Vector& fields) :
    fields_(fields)
{
}

//-----------------------------------------------------------------------------

} // namespace grid
} // namespace eckit

