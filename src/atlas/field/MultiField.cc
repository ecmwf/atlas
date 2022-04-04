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

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

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

void MultiField::initialize(const std::string& creator, const eckit::Parametrisation& params) {
    std::unique_ptr<MultiFieldCreator> MultiField_creator(MultiFieldCreatorFactory::build(creator, params));
    MultiField_creator->generate(*this, params);
}

array::Array& MultiField::allocate(array::DataType datatype, array::ArraySpec&& spec) {
    array_.reset(array::Array::create(datatype, std::move(spec)));
    return *array_;
}

//------------------------------------------------------------------------------------------------------

MultiField::MultiField() = default;

MultiField::MultiField(const std::string& creator, const eckit::Parametrisation& params) {
    initialize(creator, params);
}

const util::Metadata& MultiField::metadata() const {
    return metadata_;
}

util::Metadata& MultiField::metadata() {
    return metadata_;
}

std::vector<std::string> MultiField::field_names() const {
    std::vector<std::string> ret;
    if (fields_.size()) {
        ret.reserve(fields_.size());
    }

    for (auto it = field_index_.begin(); it != field_index_.end(); ++it) {
        ret.push_back(it->first);
    }
    return ret;
}

//-----------------------------------------------------------------------------

MultiFieldCreator::MultiFieldCreator(const eckit::Parametrisation&) {}

MultiFieldCreator::~MultiFieldCreator() = default;

MultiFieldCreator* MultiFieldCreatorFactory::build(const std::string& builder, const eckit::Parametrisation& config) {
    force_link();
    auto factory = get(builder);
    return factory->make(config);
}

//-----------------------------------------------------------------------------

}  // namespace field
}  // namespace atlas
