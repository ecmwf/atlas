/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "atlas/trans/Cache.h"

#include "pluto/pluto.h"

#include "eckit/io/DataHandle.h"

#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "atlas/runtime/Trace.h"
#include "atlas/trans/Trans.h"

namespace atlas {
namespace trans {

TransCacheFileEntry::TransCacheFileEntry(const eckit::PathName& path) {
    ATLAS_TRACE();
    Log::debug() << "Loading cache from file " << path << std::endl;
    std::unique_ptr<eckit::DataHandle> dh(path.fileHandle());

    size_ = path.size();
    data_ = pluto::host_resource()->allocate(size_, 256);

    dh->openForRead();
    dh->read(data_, size_);
    dh->close();
}

TransCacheFileEntry::~TransCacheFileEntry() {
    pluto::host_resource()->deallocate(data_, size_, 256);
}

TransCacheMemoryEntry::TransCacheMemoryEntry(const void* data, size_t size): data_(data), size_(size) {
    ATLAS_ASSERT(data_);
    ATLAS_ASSERT(size_);
}

LegendreFFTCache::LegendreFFTCache(const void* legendre_address, size_t legendre_size, const void* fft_address,
                                   size_t fft_size):
    Cache(std::make_shared<TransCacheMemoryEntry>(legendre_address, legendre_size),
          std::make_shared<TransCacheMemoryEntry>(fft_address, fft_size)) {}

LegendreFFTCache::LegendreFFTCache(const eckit::PathName& legendre_path, const eckit::PathName& fft_path):
    Cache(std::shared_ptr<TransCacheEntry>(new TransCacheFileEntry(legendre_path)),
          std::shared_ptr<TransCacheEntry>(new TransCacheFileEntry(fft_path))) {}

LegendreCache::LegendreCache(const eckit::PathName& path):
    Cache(std::shared_ptr<TransCacheEntry>(new TransCacheFileEntry(path))) {}

LegendreCache::LegendreCache(size_t size): Cache(std::make_shared<TransCacheOwnedMemoryEntry>(size)) {}

LegendreCache::LegendreCache(const void* address, size_t size):
    Cache(std::make_shared<TransCacheMemoryEntry>(address, size)) {}

Cache::Cache(const std::shared_ptr<TransCacheEntry>& legendre):
    trans_(nullptr), legendre_(legendre), fft_(new EmptyCacheEntry()) {}

Cache::Cache(const std::shared_ptr<TransCacheEntry>& legendre, const std::shared_ptr<TransCacheEntry>& fft):
    trans_(nullptr), legendre_(legendre), fft_(fft) {}

Cache::Cache(const TransImpl* trans): trans_(trans), legendre_(new EmptyCacheEntry()), fft_(new EmptyCacheEntry()) {}

Cache::Cache(): trans_(nullptr), legendre_(new EmptyCacheEntry()), fft_(new EmptyCacheEntry()) {}

Cache::operator bool() const {
    return trans_ || bool(legendre());
}

Cache::~Cache() = default;

TransCache::TransCache(const Trans& trans): Cache(trans.get()) {}

TransCacheOwnedMemoryEntry::TransCacheOwnedMemoryEntry(size_t size): size_(size) {
    if (size_) {
        data_ = pluto::host_resource()->allocate(size_, 256);
    }
}

TransCacheOwnedMemoryEntry::~TransCacheOwnedMemoryEntry() {
    if (size_) {
        pluto::host_resource()->deallocate(data_, size_, 256);
    }
}

}  // namespace trans
}  // namespace atlas
