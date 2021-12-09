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

#include "atlas/projection/Projection.h"
#include "atlas/util/Factory.h"

namespace eckit {
class Parametrisation;
}

namespace atlas {
namespace projection {

//----------------------------------------------------------------------------------------------------------------------

class ProjectionFactory : public util::Factory<ProjectionFactory> {
public:
    static std::string className() { return "ProjectionFactory"; }
    static const Projection::Implementation* build(const std::string&);
    static const Projection::Implementation* build(const std::string&, const eckit::Parametrisation&);
    using Factory::Factory;

private:
    virtual const Projection::Implementation* make(const eckit::Parametrisation&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class ProjectionBuilder : public ProjectionFactory {
private:
    virtual const Projection::Implementation* make(const eckit::Parametrisation& param) { return new T(param); }

public:
    using ProjectionFactory::ProjectionFactory;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace projection
}  // namespace atlas
