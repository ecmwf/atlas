/*
 * (C) Crown Copyright 2021 Met Office.
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

class CubedSphereTiles;

//----------------------------------------------------------------------------------------------------------------------

class CubedSphereTilesFactory : public util::Factory<CubedSphereTilesFactory> {
public:
    static std::string className() { return "CubedSphereTilesFactory"; }
    static const CubedSphereTiles* build(const std::string&);
    static const CubedSphereTiles* build(const std::string&, const eckit::Parametrisation&);
    using Factory::Factory;

private:
    virtual const CubedSphereTiles* make(const eckit::Parametrisation&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class CubedSphereTilesBuilder : public CubedSphereTilesFactory {
private:
    virtual const CubedSphereTiles* make(const eckit::Parametrisation& param) { return new T(param); }

public:
    using CubedSphereTilesFactory::CubedSphereTilesFactory;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace grid
}  // namespace atlas
