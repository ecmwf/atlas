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

#include "pluto/pluto.h"

#include "atlas/array/ArrayDataStore.h"
#include "atlas/library/Library.h"
#include "atlas/library/config.h"
#include "atlas/parallel/acc/acc.h"
#include "atlas/runtime/Exception.h"
#include "atlas/runtime/Log.h"
#include "eckit/log/Bytes.h"

#define ATLAS_ACC_DEBUG 0

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
    DataStore(size_t size): size_(size),
        host_allocator_{pluto::host_resource()},
        device_allocator_{pluto::device_resource()} {
        allocateHost();
        initialise(host_data_, size_);
        if (ATLAS_HAVE_GPU && pluto::devices()) {
            device_updated_ = false;
        }
        else {
            device_data_ = host_data_;
        }
    }

    virtual ~DataStore() {
        deallocateDevice();
        deallocateHost();
    }

    void updateDevice() const override {
        if (ATLAS_HAVE_GPU && pluto::devices()) {
            if (not device_allocated_) {
                allocateDevice();
            }
            pluto::copy_host_to_device(device_data_, host_data_, size_);
            device_updated_ = true;
        }
    }

    void updateHost() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                pluto::copy_device_to_host(host_data_, device_data_, size_);
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
        if (ATLAS_HAVE_GPU && pluto::devices()) {
            if (device_allocated_) {
                return;
            }
            if (size_) {
                device_data_ = device_allocator_.allocate(size_);
                device_allocated_ = true;
                accMap();
            }
        }
    }

    void deallocateDevice() const override {
        if (device_allocated_) {
            accUnmap();
            device_allocator_.deallocate(device_data_,size_);
            device_data_ = nullptr;
            device_allocated_ = false;
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
        if (not acc_mapped_ && acc::devices()) {
            ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
            ATLAS_ASSERT(!atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            if constexpr(ATLAS_ACC_DEBUG) {
                std::cout << "               + acc_map_data(hostptr:"<<host_data_<<", device:"<<device_data_<<", bytes:"<<footprint()<<")" <<std::endl;
            }
            atlas::acc::map((void*)host_data_, (void*)device_data_, size_ * sizeof(Value));
            ATLAS_ASSERT(atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            ATLAS_ASSERT(atlas::acc::deviceptr(host_data_) == device_data_);
            acc_mapped_ = true;
        }
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
        if (acc_mapped_) {
            if constexpr(ATLAS_ACC_DEBUG) {
                std::cout << "               - acc_unmap_data(hostptr:"<<host_data_<<", device:"<<device_data_<<", bytes:"<<footprint()<<")" <<std::endl;
            }
            ATLAS_ASSERT(atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            atlas::acc::unmap(host_data_);
            ATLAS_ASSERT(!atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            acc_mapped_ = false;
        }
    }


private:
    [[noreturn]] void throw_AllocationFailed(size_t bytes, const eckit::CodeLocation& loc) {
        std::ostringstream ss;
        ss << "AllocationFailed: Could not allocate " << eckit::Bytes(bytes);
        throw_Exception(ss.str(), loc);
    }

    void allocateHost() {
        if (size_ > 0) {
            MemoryHighWatermark::instance() += footprint();
            host_data_ = host_allocator_.allocate(size_);
         }
         else {
            host_data_ = nullptr;
         }
    }

    void deallocateHost() {
        if (host_data_) {
            host_allocator_.deallocate(host_data_, size_);
            host_data_ = nullptr;
            MemoryHighWatermark::instance() -= footprint();
        }
    }

    size_t footprint() const { return sizeof(Value) * size_; }

    size_t size_;
    Value* host_data_;
    mutable Value* device_data_{nullptr};

    mutable bool host_updated_{true};
    mutable bool device_updated_{true};
    mutable bool device_allocated_{false};
    mutable bool acc_mapped_{false};

    pluto::allocator<Value> host_allocator_;
    mutable pluto::allocator<Value> device_allocator_;
};

//------------------------------------------------------------------------------

template <typename Value>
class WrappedDataStore : public ArrayDataStore {
public:

    void init_device() {
        if (ATLAS_HAVE_GPU && pluto::devices()) {
            device_updated_ = false;
        }
        else {
            device_data_ = host_data_;
        }
    }

    WrappedDataStore(Value* host_data, size_t size): host_data_(host_data), size_(size),
        device_allocator_{pluto::device_resource()} {
        init_device();
    }

    WrappedDataStore(Value* host_data, const ArraySpec& spec):
        host_data_(host_data),
        size_(spec.size()),
        device_allocator_{pluto::device_resource()} {
        init_device();
        contiguous_ = spec.contiguous();
        if (! contiguous_) {
            int break_idx = 0;
            size_t shp_mult_rhs = spec.shape()[spec.rank() - 1];
            for (int i = spec.rank() - 1; i > 0; i--) {
                if (shp_mult_rhs != spec.strides()[i - 1]) {
                    break_idx = i - 1;
                    break;
                }
                shp_mult_rhs *= spec.shape()[i - 1];
            }
            size_t shp_mult_lhs = spec.shape()[0];
            for (int i = 1; i <= break_idx; i++) {
                shp_mult_lhs *= spec.shape()[i];
            }
            if (spec.strides()[spec.rank() - 1] > 1) {
                memcpy_h2d_pitch_ = 1;
                memcpy_d2h_pitch_ = spec.strides()[spec.rank() - 1];
                memcpy_width_ = 1;
                memcpy_height_ = spec.shape()[0] * spec.device_strides()[0];
            }
            else {
                memcpy_h2d_pitch_ = shp_mult_rhs;
                memcpy_d2h_pitch_ = spec.strides()[break_idx];
                memcpy_width_ = shp_mult_rhs;
                memcpy_height_ = shp_mult_lhs;
            }
        }
    }

    virtual ~WrappedDataStore() {
        deallocateDevice();
    }

    void updateDevice() const override {
        if (ATLAS_HAVE_GPU && pluto::devices()) {
            if (not device_allocated_) {
                allocateDevice();
            }
            if (contiguous_) {
                pluto::copy_host_to_device(device_data_, host_data_, size_);
            }
            else {
                pluto::copy_host_to_device_2D(
                    device_data_, memcpy_h2d_pitch_,
                    host_data_, memcpy_d2h_pitch_,
                    memcpy_width_, memcpy_height_);
            }
            device_updated_ = true;
        }
    }

    void updateHost() const override {
        if constexpr (ATLAS_HAVE_GPU) {
            if (device_allocated_) {
                if (contiguous_) {
                    pluto::copy_device_to_host(host_data_, device_data_, size_);
                }
                else {
                    pluto::copy_device_to_host_2D(
                        host_data_, memcpy_d2h_pitch_ ,
                        device_data_, memcpy_h2d_pitch_,
                        memcpy_width_, memcpy_height_);
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
        if (ATLAS_HAVE_GPU && pluto::devices()) {
            if (device_allocated_) {
                return;
            }
            if (size_) {
                device_data_ = device_allocator_.allocate(size_);
                device_allocated_ = true;
                if (contiguous_) {
                    accMap();
                }
            }
        }
    }

    void deallocateDevice() const override {
        if (device_allocated_) {
            if (contiguous_) {
                accUnmap();
            }
            device_allocator_.deallocate(device_data_, size_);
            device_data_ = nullptr;
            device_allocated_ = false;
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
        if (not acc_mapped_ && acc::devices()) {
            ATLAS_ASSERT(deviceAllocated(),"Could not accMap as device data is not allocated");
            ATLAS_ASSERT(!atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            if constexpr(ATLAS_ACC_DEBUG) {
                std::cout << "               + acc_map_data(hostptr:"<<host_data_<<", device:"<<device_data_<<", bytes:"<<size_ * sizeof(Value)<<")" <<std::endl;
            }
            atlas::acc::map((void*)host_data_, (void*)device_data_, size_ * sizeof(Value));
            ATLAS_ASSERT(atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            ATLAS_ASSERT(atlas::acc::deviceptr(host_data_) == device_data_);
            acc_mapped_ = true;
        }
    }

    bool accMapped() const override {
        return acc_mapped_;
    }

    void accUnmap() const override {
        if (acc_mapped_) {
            ATLAS_ASSERT(atlas::acc::is_present(host_data_, size_ * sizeof(Value)));
            if constexpr(ATLAS_ACC_DEBUG) {
                std::cout << "               - acc_unmap_data(hostptr:"<<host_data_<<", device:"<<device_data_<<", bytes:"<<size_ * sizeof(Value)<<")" <<std::endl;
            }
            atlas::acc::unmap(host_data_);
            acc_mapped_ = false;
        }
    }

private:
    Value* host_data_;
    size_t size_;
    mutable Value* device_data_;

    bool contiguous_{true};
    size_t memcpy_h2d_pitch_;
    size_t memcpy_d2h_pitch_;
    size_t memcpy_height_;
    size_t memcpy_width_;

    mutable bool host_updated_{true};
    mutable bool device_updated_{true};
    mutable bool device_allocated_{false};
    mutable bool acc_mapped_{false};

    mutable pluto::allocator<Value> device_allocator_;
};

}  // namespace native
}  // namespace array
}  // namespace atlas
