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
namespace spacing {

//----------------------------------------------------------------------------------------------------------------------

class Spacing;
class SpacingFactory : public util::Factory<SpacingFactory> {
public:
    static std::string className() { return "SpacingFactory"; }
    static const Spacing* build(const std::string&);
    static const Spacing* build(const std::string&, const eckit::Parametrisation&);
    using Factory::Factory;

private:
    virtual const Spacing* make(const eckit::Parametrisation&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class SpacingBuilder : public SpacingFactory {
private:
    virtual const Spacing* make(const eckit::Parametrisation& param) { return new T(param); }

public:
    using SpacingFactory::SpacingFactory;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace spacing
}  // namespace grid
}  // namespace atlas
