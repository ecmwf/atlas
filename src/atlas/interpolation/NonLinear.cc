/*
 * (C) Copyright 1996- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#include "atlas/interpolation/NonLinear.h"

#include "atlas/interpolation/nonlinear/NonLinear.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Config.h"


namespace atlas {
namespace interpolation {


namespace {
std::string config_type(const eckit::Parametrisation& config) {
    std::string value;
    ATLAS_ASSERT_MSG(config.get("type", value), "");
    return value;
}
}  // namespace


NonLinear::NonLinear(): Handle(nullptr) {}


NonLinear::NonLinear(const NonLinear::Config& config):
    Handle(nonlinear::NonLinearFactory::build(config_type(config), config)) {}


NonLinear::NonLinear(const std::string& type, const NonLinear::Config& config):
    Handle(nonlinear::NonLinearFactory::build(type, config)) {}


bool NonLinear::operator()(const Field& f) const {
    return operator bool() && get()->applicable(f);
}


bool NonLinear::execute(NonLinear::Matrix& W, const Field& f) const {
    ATLAS_ASSERT_MSG(operator bool(), "NonLinear: ObjectHandle not setup");
    return get()->execute(W, f);
}

bool NonLinear::execute(NonLinear::Matrix& W, const array::Array& a, const NonLinear::Config& c) const {
    ATLAS_ASSERT_MSG(operator bool(), "NonLinear: ObjectHandle not setup");
    return get()->execute(W, a, c);
}

}  // namespace interpolation
}  // namespace atlas
