/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "TraceMemoryResource.h"

#include <limits>

namespace pluto {

int TraceMemoryResource::nest = 0;

void* TraceMemoryResource::do_allocate(std::size_t bytes, std::size_t alignment) {
    nest++;
    if (TraceOptions::instance().enabled) {
        *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (alloc)] " << bytes << " bytes with alignment " << alignment << std::endl;
    }
    auto ptr = mr_->allocate(bytes, alignment);
    nest--;
    return  ptr;
}

void TraceMemoryResource::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    nest++;
    // print_stacktrace(out_, 4*nest, CURRENT_STACKTRACE());
    if (TraceOptions::instance().enabled) {
        if (bytes == std::numeric_limits<std::size_t>::max()) {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc)] " << p << std::endl;
        }
        else {
            *TraceOptions::instance().out << std::string(4*nest,' ') << "[" << name_ << " (dealloc)] " << p << " : " << bytes << " bytes with alignment " << alignment << std::endl;
        }
    }
    mr_->deallocate(p, bytes, alignment);
    nest--;
}

}
