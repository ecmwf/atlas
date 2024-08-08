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

#include "hic/hic.h"


#if ATLAS_HAVE_ACC
#include "atlas_acc_support/atlas_acc_map_data.h"
#endif

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
        allocateHost();
        initialise(host_data_, size_);
        if constexpr (ATLAS_HAVE_GPU) {
            device_updated_ = false;
        }
        else {
            device_data_ = host_data_;
        }
    }

    ~DataStore() {
        deallocateDevice();
        deallocateHost();
    }

    void updateDevice() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (not device_allocated_) {
                allocateDevice();
            }
            hicError_t err = hicMemcpy(device_data_, host_data_, size_*sizeof(Value), hicMemcpyHostToDevice);
            if (err != hicSuccess) {
                throw_AssertionFailed("Failed to updateDevice: "+std::string(hicGetErrorString(err)), Here());
            }
            device_updated_ = true;
        }
    }

    void updateHost() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                hicError_t err = hicMemcpy(host_data_, device_data_, size_*sizeof(Value), hicMemcpyDeviceToHost);
                if (err != hicSuccess) {
                    throw_AssertionFailed("Failed to updateHost: "+std::string(hicGetErrorString(err)), Here());
                }
                host_updated_ = true;
            }
        }
    }

    bool valid() const override { return true; }

    void syncHostDevice() const override {
        if (host_updated_ and device_updated_) {
            return; // nothing to do
        }
        if (not (host_updated_ or device_updated_)) {
            throw_AssertionFailed("syncHostDevice() could not figure out which of host or device is up to date. "
                                  "Probably it was forgotten to use setDeviceNeedsUpdate(true) or setDeviceNeedsUpdate(true)",
                                  Here());
        }

        if (not device_updated_) {
            updateDevice();
        }
        else if (not host_updated_) {
            updateHost();
        }
    }

    bool deviceAllocated() const override { return device_allocated_; }

    void allocateDevice() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                return;
            }
            if (size_) {
                hicError_t err = hicMalloc((void**)&device_data_, sizeof(Value)*size_);
                if (err != hicSuccess) {
                    throw_AssertionFailed("Failed to allocate GPU memory: " + std::string(hicGetErrorString(err)), Here());
                }
                device_allocated_ = true;
                accMap();
            }
        }
    }

    void deallocateDevice() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                accUnmap();
                hicError_t err = hicFree(device_data_);
                if (err != hicSuccess) {
                    throw_AssertionFailed("Failed to deallocate GPU memory: " + std::string(hicGetErrorString(err)), Here());
                }
                device_data_ = nullptr;
                device_allocated_ = false;
            }
        }
    }

    bool hostNeedsUpdate() const override { return (not host_updated_); }

    bool deviceNeedsUpdate() const override { return (not device_updated_); }

    void setHostNeedsUpdate(bool v) const override { host_updated_ = (not v); }

    void setDeviceNeedsUpdate(bool v) const override { device_updated_ = (not v); }

    void reactivateDeviceWriteViews() const override {}

    void reactivateHostWriteViews() const override {}

    void* voidDataStore() override { return static_cast<void*>(host_data_); }

    void* voidHostData() override { return static_cast<void*>(host_data_); }

    void* voidDeviceData() override { return static_cast<void*>(device_data_); }

    void accMap() const override {
#if ATLAS_HAVE_ACC
        if (not acc_mapped_) {
            ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
            atlas_acc_map_data((void*)host_data_, (void*)device_data_, size_ * sizeof(Value));
            acc_mapped_ = true;
        }
#endif
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
#if ATLAS_HAVE_ACC
        if (acc_mapped_) {
            atlas_acc_unmap_data(host_data_);
            acc_mapped_ = false;
        }
#endif
    }


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
        if (ptr) {
            free(ptr);
            ptr = nullptr;
            MemoryHighWatermark::instance() -= footprint();
        }
    }

    void allocateHost() {
        alloc_aligned(host_data_, size_);
    }

    void deallocateHost() {
        free_aligned(host_data_);
    }

    size_t footprint() const { return sizeof(Value) * size_; }

    size_t size_;
    Value* host_data_;
    mutable Value* device_data_{nullptr};

    mutable bool host_updated_{true};
    mutable bool device_updated_{true};
    mutable bool device_allocated_{false};
    mutable bool acc_mapped_{false};

};

//------------------------------------------------------------------------------

template <typename Value>
class WrappedDataStore : public ArrayDataStore {
public:
    WrappedDataStore(Value* host_data, size_t size): host_data_(host_data), size_(size) {
        if constexpr (ATLAS_HAVE_GPU) {
            device_updated_ = false;
        }
        else {
            device_data_ = host_data_;
        }
    }

    void updateDevice() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (not device_allocated_) {
                allocateDevice();
            }
            hicError_t err = hicMemcpy(device_data_, host_data_, size_*sizeof(Value), hicMemcpyHostToDevice);
            if (err != hicSuccess) {
                throw_AssertionFailed("Failed to updateDevice: "+std::string(hicGetErrorString(err)), Here());
            }
            device_updated_ = true;
        }
    }

    void updateHost() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                hicError_t err = hicMemcpy(host_data_, device_data_, size_*sizeof(Value), hicMemcpyDeviceToHost);
                if (err != hicSuccess) {
                    throw_AssertionFailed("Failed to updateHost: "+std::string(hicGetErrorString(err)), Here());
                }
                host_updated_ = true;
            }
        }
    }

    bool valid() const override { return true; }

    void syncHostDevice() const override {
        if (host_updated_ and device_updated_) {
            return; // nothing to do
        }
        if (not (host_updated_ or device_updated_)) {
            throw_AssertionFailed("syncHostDevice() could not figure out which of host or device is up to date. "
                                  "Probably it was forgotten to use setDeviceNeedsUpdate(true) or setDeviceNeedsUpdate(true)",
                                  Here());
        }

        if (not device_updated_) {
            updateDevice();
        }
        else if (not host_updated_) {
            updateHost();
        }
    }

    bool deviceAllocated() const override { return device_allocated_; }

    void allocateDevice() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                return;
            }
            if (size_) {
                hicError_t err = hicMalloc((void**)&device_data_, sizeof(Value)*size_);
                if (err != hicSuccess) {
                    throw_AssertionFailed("Failed to allocate GPU memory: " + std::string(hicGetErrorString(err)), Here());
                }
                device_allocated_ = true;
                accMap();
            }
        }
    }

    void deallocateDevice() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                accUnmap();
                hicError_t err = hicFree(device_data_);
                if (err != hicSuccess) {
                    throw_AssertionFailed("Failed to deallocate GPU memory: " + std::string(hicGetErrorString(err)), Here());
                }
                device_data_ = nullptr;
                device_allocated_ = false;
            }
        }
    }

    bool hostNeedsUpdate() const override { return (not host_updated_); }

    bool deviceNeedsUpdate() const override { return (not device_updated_); }

    void setHostNeedsUpdate(bool v) const override { host_updated_ = (not v); }

    void setDeviceNeedsUpdate(bool v) const override { device_updated_ = (not v); }

    void reactivateDeviceWriteViews() const override {}

    void reactivateHostWriteViews() const override {}

    void* voidDataStore() override { return static_cast<void*>(host_data_); }

    void* voidHostData() override { return static_cast<void*>(host_data_); }

    void* voidDeviceData() override { return static_cast<void*>(device_data_); }

    void accMap() const override {
#if ATLAS_HAVE_ACC
        if (not acc_mapped_) {
            ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
            atlas_acc_map_data((void*)host_data_, (void*)device_data_, size_ * sizeof(Value));
            acc_mapped_ = true;
        }
#endif
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
#if ATLAS_HAVE_ACC
        if (acc_mapped_) {
            atlas_acc_unmap_data(host_data_);
            acc_mapped_ = false;
        }
#endif
    }

private:
    size_t size_;
    Value* host_data_;
    mutable Value* device_data_;

    mutable bool host_updated_{true};
    mutable bool device_updated_{true};
    mutable bool device_allocated_{false};
    mutable bool acc_mapped_{false};
};

}  // namespace native
}  // namespace array
}  // namespace atlas
