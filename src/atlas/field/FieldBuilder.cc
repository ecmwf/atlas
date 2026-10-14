/*
 * (C) Copyright 2026 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/field/FieldBuilder.h"

#include "atlas/field/builders/GridBoxArea.h"

namespace atlas {
namespace field {

namespace {

void force_link() {
    static struct Link {
        Link() { FieldBuilderFactoryBuilder<builders::GridBoxArea>(); }
    } link;
}

}  // namespace

const AbstractFieldBuilder* FieldBuilderFactory::build(const util::Config& config,
                                                       const util::ObjectHandleBase& object) {
    force_link();
    const std::string type = config.getString("type");
    return get(type)->make(config, object);
}

}  // namespace field
}  // namespace atlas