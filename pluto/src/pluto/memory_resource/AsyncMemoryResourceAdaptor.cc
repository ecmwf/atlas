/*
 * (C) Copyright 2024- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "AsyncMemoryResourceAdaptor.h"

#include "pluto/stream.h"

namespace pluto {

void* AsyncMemoryResourceAdaptor::do_allocate(std::size_t bytes, std::size_t alignment) {
    if (async_mr_) {
        return async_mr_->allocate_async(bytes, alignment, get_current_stream());
    }
    else {
        return mr_->allocate(bytes, alignment);
    }
}

void AsyncMemoryResourceAdaptor::do_deallocate(void* p, std::size_t bytes, std::size_t alignment) {
    if (async_mr_) {
        async_mr_->deallocate_async(p, bytes, alignment, get_current_stream());
    }
    else {
        mr_->deallocate(p, bytes, alignment);
    }
}

}  // namespace pluto
