/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

// file deepcode ignore CppMemoryLeak: static pointers for global registry are OK and will be cleaned up at end

#include "atlas/field/MultiFieldCreator.h"

#include <map>
#include <sstream>

#include "eckit/thread/AutoLock.h"
#include "eckit/thread/Mutex.h"

#include "atlas/field/MultiFieldCreatorIFS.h"
#include "atlas/field/MultiFieldCreatorArray.h"
#include "atlas/grid/Grid.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"


namespace atlas {
namespace field {


namespace {

void force_link() {
  static struct Link {
    Link() {
      MultiFieldCreatorBuilder<MultiFieldCreatorIFS>();
      MultiFieldCreatorBuilder<MultiFieldCreatorArray>();
    }
  } link;
}

}  // namespace

// ------------------------------------------------------------------

MultiFieldCreator::MultiFieldCreator() = default;

MultiFieldCreator::~MultiFieldCreator() = default;

MultiFieldCreator* MultiFieldCreatorFactory::build(const std::string& builder, const eckit::Configuration& config) {
    force_link();
    auto factory = get(builder);
    return factory->make(config);
}

}  // namespace field
}  // namespace atlas
