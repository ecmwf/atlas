/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "MultiField.h"

#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <mutex>

#include "atlas/field/Field.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"

namespace atlas {
namespace field {

namespace {
void force_link() {
    static struct Link {
        Link() {
            ;
            // For static linking add here something like
            // MultiFieldCreatorBuilder<T>();
        }
    } link;
}
}  // namespace

//-----------------------------------------------------------------------------

class MultiFieldArrayRegistry : public field::FieldObserver {
private:
    MultiFieldArrayRegistry() {}

public:
    static MultiFieldArrayRegistry& instance() {
        static MultiFieldArrayRegistry inst;
        return inst;
    }
    void onFieldDestruction(FieldImpl& field) override {
        std::lock_guard<std::mutex> guard(lock_);
        map_.erase(&field);
    }

    ~MultiFieldArrayRegistry() override = default;

    void add(Field& field, std::shared_ptr<array::Array> array) {
        std::lock_guard<std::mutex> guard(lock_);
        map_.emplace(field.get(),array);
        field->attachObserver(*this);
    }

public:
    std::mutex lock_;
    std::map<FieldImpl*,std::shared_ptr<array::Array>> map_;

};

MultiFieldCreator::MultiFieldCreator(const eckit::Configuration&) {}

MultiFieldCreator::~MultiFieldCreator() = default;

MultiFieldCreator* MultiFieldCreatorFactory::build(const std::string& builder, const eckit::Configuration& config) {
    force_link();
    auto factory = get(builder);
    return factory->make(config);
}

MultiField::MultiField(const eckit::Configuration& config) {
    std::string type;
    if (!config.get("type", type)) {
        ATLAS_THROW_EXCEPTION("Could not find \"type\" in configuration");
    }
    std::unique_ptr<MultiFieldCreator> creator(MultiFieldCreatorFactory::build(type, config));
    reset(creator->create(config));
}

void MultiFieldImpl::add(Field& field) {
    ATLAS_ASSERT(not fieldset_.has(field.name()), "Field with name \"" + field.name() + "\" already exists!");
    fieldset_.add(field);
    MultiFieldArrayRegistry::instance().add(field,array_);

}

//-----------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
