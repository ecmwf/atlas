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

#include <iosfwd>
#include <string>
#include <vector>

#include "atlas/interpolation/method/Method.h"
#include "atlas/util/Factory.h"
#include "atlas/util/Object.h"

#include "eckit/config/Configuration.h"

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
}  // namespace atlas

namespace atlas {
namespace interpolation {

struct MethodFactory : public util::Factory<MethodFactory> {
public:
    static std::string className() { return "MethodFactory"; }
    static Method* build(const std::string& name, const Method::Config&);

protected:
    virtual Method* make(const Method::Config&) = 0;
    using Factory::Factory;
};

template <class T>
struct MethodBuilder : public MethodFactory {
    using MethodFactory::MethodFactory;

private:
    virtual Method* make(const Method::Config& config) { return new T(config); }
};

}  // namespace interpolation
}  // namespace atlas
