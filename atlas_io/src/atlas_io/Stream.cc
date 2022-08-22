/*
 * (C) Copyright 2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas_io/Stream.h"

#include "eckit/io/DataHandle.h"

#include "atlas_io/atlas_compat.h"
#include "atlas_io/detail/Assert.h"

namespace atlas {
namespace io {

//---------------------------------------------------------------------------------------------------------------------

Stream::Stream(eckit::DataHandle& datahandle): ptr_(&datahandle) {}

Stream::Stream(eckit::DataHandle* datahandle): shared_(datahandle), ptr_(shared_.get()) {}

Stream::Stream(std::shared_ptr<eckit::DataHandle> datahandle): shared_(datahandle), ptr_(shared_.get()) {}

Stream::Stream(const Stream& other) = default;

eckit::DataHandle& Stream::datahandle() {
    ATLAS_IO_ASSERT(ptr_ != nullptr);
    return *ptr_;
}

uint64_t Stream::seek(uint64_t offset) {
    ATLAS_IO_ASSERT(ptr_ != nullptr);
    return std::uint64_t(ptr_->seek(static_cast<long long>(offset)));
}

uint64_t Stream::position() {
    ATLAS_IO_ASSERT(ptr_ != nullptr);
    return std::uint64_t(ptr_->position());
}

uint64_t Stream::write(const void* data, size_t length) {
    ATLAS_IO_ASSERT(ptr_ != nullptr);
    return std::uint64_t(ptr_->write(data, static_cast<long>(length)));
}

uint64_t Stream::read(void* data, size_t length) {
    ATLAS_IO_ASSERT(ptr_ != nullptr);
    return std::uint64_t(ptr_->read(data, static_cast<long>(length)));
}

//---------------------------------------------------------------------------------------------------------------------

}  // namespace io
}  // namespace atlas
