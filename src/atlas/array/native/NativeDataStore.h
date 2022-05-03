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

#include <algorithm>  // std::fill
#include <atomic>
#include <cstdlib>  // posix_memalign
#include <limits>   // std::numeric_limits<T>::signaling_NaN
#include <sstream>

#include "atlas/array/ArrayDataStore.h"
#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "eckit/log/Bytes.h"

//------------------------------------------------------------------------------

namespace atlas {
namespace array {
namespace native {

struct MemoryHighWatermark {
    std::atomic<size_t> bytes_{0};
    std::atomic<size_t> high_{0};
    void print(std::ostream& out) const { out << eckit::Bytes(double(bytes_)); }
    friend std::ostream& operator<<(std::ostream& out, const MemoryHighWatermark& v) {
        v.print(out);
        return out;
    }
    MemoryHighWatermark& operator+=(const size_t& bytes) {
        bytes_ += bytes;
        update_maximum();
        if (atlas::Library::instance().traceMemory()) {
            Log::trace() << "Memory: " << eckit::Bytes(double(bytes_)) << "\t( +" << eckit::Bytes(double(bytes))
                         << " \t| high watermark " << eckit::Bytes(double(high_)) << "\t)" << std::endl;
        }
        return *this;
    }
    MemoryHighWatermark& operator-=(const size_t& bytes) {
        bytes_ -= bytes;
        if (atlas::Library::instance().traceMemory()) {
            Log::trace() << "Memory: " << eckit::Bytes(double(bytes_)) << "\t( -" << eckit::Bytes(double(bytes))
                         << " \t| high watermark " << eckit::Bytes(double(high_)) << "\t)" << std::endl;
        }
        return *this;
    }

private:
    MemoryHighWatermark() = default;
    void update_maximum() noexcept {
        size_t prev_value = high_;
        while (prev_value < bytes_ && !high_.compare_exchange_weak(prev_value, bytes_)) {
        }
    }

public:
    static MemoryHighWatermark& instance() {
        static MemoryHighWatermark _instance;
        return _instance;
    }
};

template <typename Value>
static constexpr Value invalid_value() {
    return std::numeric_limits<Value>::has_signaling_NaN ? std::numeric_limits<Value>::signaling_NaN()
           : std::numeric_limits<Value>::has_quiet_NaN   ? std::numeric_limits<Value>::quiet_NaN()
           : std::numeric_limits<Value>::has_infinity    ? std::numeric_limits<Value>::infinity()
                                                         : std::numeric_limits<Value>::max();
}

#if ATLAS_INIT_SNAN
template <typename Value>
void initialise(Value array[], size_t size) {
    std::fill_n(array, size, invalid_value<Value>());
}
#else
template <typename Value>
void initialise(Value[], size_t) {}
#endif

template <typename Value>
class DataStore : public ArrayDataStore {
public:
    DataStore(size_t size): size_(size) {
        alloc_aligned(data_store_, size_);
        initialise(data_store_, size_);
    }

    virtual ~DataStore() override { free_aligned(data_store_); }

    virtual void updateDevice() const override {}

    virtual void updateHost() const override {}

    virtual bool valid() const override { return true; }

    virtual void syncHostDevice() const override {}

    virtual bool deviceAllocated() const override { return false; }

    virtual void allocateDevice() const override {}

    virtual void deallocateDevice() const override {}

    virtual bool hostNeedsUpdate() const override { return false; }

    virtual bool deviceNeedsUpdate() const override { return false; }

    virtual void setHostNeedsUpdate(bool) const override {}

    virtual void setDeviceNeedsUpdate(bool) const override {}

    virtual void* voidDataStore() override { return static_cast<void*>(data_store_); }

    virtual void* voidHostData() override { return static_cast<void*>(data_store_); }

    virtual void* voidDeviceData() override { return static_cast<void*>(data_store_); }

private:
    [[noreturn]] void throw_AllocationFailed(size_t bytes, const eckit::CodeLocation& loc) {
        std::ostringstream ss;
        ss << "AllocationFailed: Could not allocate " << eckit::Bytes(bytes);
        throw_Exception(ss.str(), loc);
    }

    void alloc_aligned(Value*& ptr, size_t n) {
        if (n > 0) {
            const size_t alignment = 64 * sizeof(Value);
            size_t bytes           = sizeof(Value) * n;
            MemoryHighWatermark::instance() += bytes;

            int err = posix_memalign((void**)&ptr, alignment, bytes);
            if (err) {
                throw_AllocationFailed(bytes, Here());
            }
        }
        else {
            ptr = nullptr;
        }
    }

    void free_aligned(Value*& ptr) {
        if (size_) {
            free(ptr);
            ptr = nullptr;
            MemoryHighWatermark::instance() -= footprint();
        }
    }

    size_t footprint() const { return sizeof(Value) * size_; }

    Value* data_store_;
    size_t size_;
};

//------------------------------------------------------------------------------

template <typename Value>
class WrappedDataStore : public ArrayDataStore {
public:
    WrappedDataStore(Value* data_store): data_store_(data_store) {}

    virtual void updateHost() const override {}

    virtual void updateDevice() const override {}

    virtual bool valid() const override { return true; }

    virtual void syncHostDevice() const override {}

    virtual bool deviceAllocated() const override { return false; }

    virtual void allocateDevice() const override {}

    virtual void deallocateDevice() const override {}

    virtual bool hostNeedsUpdate() const override { return true; }

    virtual bool deviceNeedsUpdate() const override { return false; }

    virtual void setHostNeedsUpdate(bool) const override {}

    virtual void setDeviceNeedsUpdate(bool) const override {}

    virtual void* voidDataStore() override { return static_cast<void*>(data_store_); }

    virtual void* voidHostData() override { return static_cast<void*>(data_store_); }

    virtual void* voidDeviceData() override { return static_cast<void*>(data_store_); }

private:
    Value* data_store_;
};

}  // namespace native
}  // namespace array
}  // namespace atlas
