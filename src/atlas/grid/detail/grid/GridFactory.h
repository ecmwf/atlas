/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"

namespace atlas {
namespace grid {

namespace detail {
namespace grid {
class Grid;
}
}  // namespace detail
using GridImpl = detail::grid::Grid;

//----------------------------------------------------------------------------------------------------------------------

class GridFactory : public util::Factory<GridFactory> {
public:
    static std::string className() { return "GridFactory"; }
    static const GridImpl* build(const std::string&, const util::Config&);
    using Factory::Factory;

private:
    virtual const GridImpl* make(const util::Config&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class GridFactoryBuilder : public GridFactory {
private:
    virtual const GridImpl* make(const util::Config& config) override { return T::create(config); }

public:
    using GridFactory::GridFactory;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace grid
}  // namespace atlas
