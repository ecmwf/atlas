/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <memory>

#include "eckit/filesystem/PathName.h"
#include "eckit/io/Buffer.h"

#include "atlas/util/ObjectHandle.h"

//-----------------------------------------------------------------------------
// Forward declarations

namespace atlas {
class Field;
class FieldSet;
class FunctionSpace;
class Grid;
class Domain;
namespace trans {
class TransImpl;
class Trans;
}  // namespace trans
}  // namespace atlas

//-----------------------------------------------------------------------------

namespace atlas {
namespace trans {

//-----------------------------------------------------------------------------

class TransCacheEntry {
public:
    operator bool() const { return size() != 0; }
    virtual ~TransCacheEntry()       = default;
    virtual size_t size() const      = 0;
    virtual const void* data() const = 0;
};

//-----------------------------------------------------------------------------

class EmptyCacheEntry final : public TransCacheEntry {
public:
    virtual size_t size() const override { return 0; }
    virtual const void* data() const override { return nullptr; }
};

//-----------------------------------------------------------------------------
class TransCacheFileEntry final : public TransCacheEntry {
private:
    eckit::Buffer buffer_;

public:
    TransCacheFileEntry(const eckit::PathName& path);
    virtual size_t size() const override { return buffer_.size(); }
    virtual const void* data() const override { return buffer_.data(); }
};

//-----------------------------------------------------------------------------

class TransCacheMemoryEntry final : public TransCacheEntry {
public:
    TransCacheMemoryEntry(const void* data, size_t size);
    virtual const void* data() const override { return data_; }
    virtual size_t size() const override { return size_; }

private:
    const void* data_;
    const size_t size_;
};

//-----------------------------------------------------------------------------

class TransCacheOwnedMemoryEntry final : public TransCacheEntry {
public:
    TransCacheOwnedMemoryEntry(size_t size);
    virtual ~TransCacheOwnedMemoryEntry() override;
    virtual const void* data() const override { return data_; }
    virtual size_t size() const override { return size_; }

private:
    void* data_        = nullptr;
    const size_t size_ = 0;
};

//-----------------------------------------------------------------------------

class Cache {
public:
    Cache();
    Cache(const Cache& other) = default;
    operator bool() const;
    const TransImpl* trans() const { return trans_.get(); }
    const TransCacheEntry& legendre() const { return *legendre_; }
    const TransCacheEntry& fft() const { return *fft_; }
    virtual ~Cache();

protected:
    Cache(const std::shared_ptr<TransCacheEntry>& legendre);
    Cache(const std::shared_ptr<TransCacheEntry>& legendre, const std::shared_ptr<TransCacheEntry>& fft);
    Cache(const TransImpl*);

private:
    util::ObjectHandle<const TransImpl> trans_;
    std::shared_ptr<TransCacheEntry> legendre_;
    std::shared_ptr<TransCacheEntry> fft_;
};

class TransCache : public Cache {
public:
    TransCache(const Trans&);
};


class LegendreCache : public Cache {
public:
    LegendreCache(size_t size);
    LegendreCache(const void* address, size_t size);
    LegendreCache(const eckit::PathName& path);
};

class LegendreFFTCache : public Cache {
public:
    LegendreFFTCache(const void* legendre_address, size_t legendre_size, const void* fft_address, size_t fft_size);
    LegendreFFTCache(const eckit::PathName& legendre_path, const eckit::PathName& fft_path);
};

//----------------------------------------------------------------------------------------------------------------------

}  // namespace trans
}  // namespace atlas
