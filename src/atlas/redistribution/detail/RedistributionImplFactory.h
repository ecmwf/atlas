/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#pragma once

#include <string>

#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"

namespace atlas {
namespace redistribution {
namespace detail {

class RedistributionImpl;

//----------------------------------------------------------------------------------------------------------------------

class RedistributionImplFactory : public util::Factory<RedistributionImplFactory> {
public:
    static std::string className() { return "RedistributionFactory"; }
    static RedistributionImpl* build(const std::string&);
    using Factory::Factory;

private:
    virtual RedistributionImpl* make() = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class RedistributionImplBuilder : public RedistributionImplFactory {
private:
    virtual RedistributionImpl* make() { return new T(); }

public:
    using RedistributionImplFactory::RedistributionImplFactory;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace redistribution
}  // namespace atlas
