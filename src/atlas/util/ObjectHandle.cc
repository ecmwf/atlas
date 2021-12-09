/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <ostream>

#include "atlas/runtime/Log.h"
#include "atlas/util/Object.h"
#include "atlas/util/ObjectHandle.h"

namespace atlas {
namespace util {

//----------------------------------------------------------------------------------------------------------------------

int ObjectHandleBase::owners() const {
    return static_cast<int>(object_->owners());
}

void ObjectHandleBase::release() {
    object_->lock();
    object_->detach(); /* lock/unlock in detach() isn't sufficient, else there is race
                               condition on owners() */

    if (object_->owners() == 0) {
        object_->unlock();
        delete object_;
        object_ = nullptr;
        return;
    }
    object_->unlock();
    object_ = nullptr;
}

void ObjectHandleBase::assign(const ObjectHandleBase& other) {
    assign(other.object_);
}

void ObjectHandleBase::assign(const Object* other) {
    if (object_) {
        release();
    }

    object_ = const_cast<Object*>(other);

    attach();
}


void ObjectHandleBase::attach() {
    if (!null()) {
        object_->attach();
    }
}

// ------------------------------------------------------------------

}  // namespace util
}  // namespace atlas
