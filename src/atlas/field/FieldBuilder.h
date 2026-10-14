/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <string>

#include "atlas/field/Field.h"
#include "atlas/util/Config.h"
#include "atlas/util/Factory.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
namespace field {

/// Polymorphic implementation interface for operators that construct a Field.
/// Implementations retain any source object state needed by the call operator.
class AbstractFieldBuilder : public util::Object {
public:
    virtual Field operator()() const = 0;
};

/// Registry of named field-builder implementations.
/// The configuration `type` selects a factory, which receives a generic Atlas object handle and the complete configuration.
class FieldBuilderFactory : public util::Factory<FieldBuilderFactory> {
public:
    static std::string className() { return "FieldBuilderFactory"; }

    static const AbstractFieldBuilder* build(const util::Config& config, const util::ObjectHandleBase& object);

    using Factory::Factory;

private:
    virtual const AbstractFieldBuilder* make(const util::Config&, const util::ObjectHandleBase&) = 0;
};

/// Factory adapter that self-registers a concrete field-builder implementation.
/// ConcreteFieldBuilder must accept `(const util::Config&, const util::ObjectHandleBase&)`.
template <typename ConcreteFieldBuilder>
class FieldBuilderFactoryBuilder : public FieldBuilderFactory {
public:
    using FieldBuilderFactory::FieldBuilderFactory;

private:
    const AbstractFieldBuilder* make(const util::Config& config, const util::ObjectHandleBase& object) override {
        return new ConcreteFieldBuilder(config, object);
    }
};

/// Owning public handle for a configured field-building operator.
/// Construction resolves the implementation through FieldBuilderFactory; invocation returns the newly built Field.
class FieldBuilder : public util::ObjectHandle<AbstractFieldBuilder> {
public:
    using Handle::Handle;

    FieldBuilder(const util::Config& config, const util::ObjectHandleBase& object):
        Handle(FieldBuilderFactory::build(config, object)) {}
    FieldBuilder(const std::string& type, const util::ObjectHandleBase& object,
                 const util::Config& config = util::NoConfig()):
        Handle(FieldBuilderFactory::build(util::Config("type", type) | config, object)) {}

    Field operator()() const { return get()->operator()(); }
};

}  // namespace field
}  // namespace atlas