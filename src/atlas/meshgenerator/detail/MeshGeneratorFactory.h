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
namespace meshgenerator {

//----------------------------------------------------------------------------------------------------------------------

class MeshGeneratorImpl;
class MeshGeneratorFactory : public util::Factory<MeshGeneratorFactory> {
public:
    static std::string className() { return "MeshGeneratorFactory"; }
    static const MeshGeneratorImpl* build(const std::string&);
    static const MeshGeneratorImpl* build(const std::string&, const eckit::Parametrisation&);
    using Factory::Factory;

private:
    virtual const MeshGeneratorImpl* make(const eckit::Parametrisation&) = 0;
};

//----------------------------------------------------------------------------------------------------------------------

template <class T>
class MeshGeneratorBuilder : public MeshGeneratorFactory {
private:
    virtual const MeshGeneratorImpl* make(const eckit::Parametrisation& param) { return new T(param); }

public:
    using MeshGeneratorFactory::MeshGeneratorFactory;
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace meshgenerator
}  // namespace atlas
