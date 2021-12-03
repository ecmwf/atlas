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


#pragma once

#include <string>

#include "eckit/config/Parametrisation.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/Factory.h"
#include "atlas/util/Object.h"


namespace atlas {
class Field;
}


namespace atlas {
namespace field {
namespace detail {


/**
 * @brief Missing values indicator base class
 */
struct MissingValue : util::Object {
    using Config                        = eckit::Parametrisation;
    virtual ~MissingValue()             = default;
    virtual bool isnan() const          = 0;
    virtual void metadata(Field&) const = 0;

    virtual bool operator()(const double&) const { ATLAS_NOTIMPLEMENTED; }
    virtual bool operator()(const float&) const { ATLAS_NOTIMPLEMENTED; }
    virtual bool operator()(const int&) const { ATLAS_NOTIMPLEMENTED; }
    virtual bool operator()(const long&) const { ATLAS_NOTIMPLEMENTED; }
    virtual bool operator()(const unsigned long&) const { ATLAS_NOTIMPLEMENTED; }
};


/**
 * @brief Missing values indicator factory
 */
struct MissingValueFactory : util::Factory<MissingValueFactory> {
    using Config = MissingValue::Config;
    using Factory::Factory;
    static std::string className() { return "MissingValueFactory"; }
    static const MissingValue* build(const std::string&, const Config&);

private:
    virtual const MissingValue* make(const Config&) = 0;
};


/**
 * @brief Missing values indicator builder for factory registration
 */
template <class T>
class MissingValueFactoryBuilder : public MissingValueFactory {
private:
    virtual const MissingValue* make(const Config& config) override { return new T(config); }

public:
    using MissingValueFactory::MissingValueFactory;
};


}  // namespace detail
}  // namespace field
}  // namespace atlas
